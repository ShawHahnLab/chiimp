# Misc utility functions and variables used by the others.

# Sample data frame helpers -----------------------------------------------

#' Check Sample Data for Potential Allele Matches
#'
#' Check the entries in a processed sample data frame for potential matches to a
#' given locus.
#'
#' @param sample_data data frame of processed data for sample as produced by
#'   [analyze_sample].
#' @param locus_name character name of locus to match against.
#'
#' @return logical vector of entries for potential alleles.
#' @md
full_locus_match <- function(sample_data, locus_name) {
  with(sample_data,
       as.character(MatchingLocus) == locus_name &
         MotifMatch &
         LengthMatch)
}

#' Create identifiers for STR Data
#'
#' Create entry IDs for the given data frame, using whichever STR-related
#' metadata columns are available.  These will not necessarily be unique.
#'
#' @param data STR data frame, such as produced by [prepare_dataset] or a
#'   summary from [summarize_dataset].
#'
#' @return character vector of entry identifiers
#' @md
make_entry_id <- function(data) {
  cols_names <- c("Dataset", "Sample", "Name", "Replicate", "Locus")
  cols_idx <- match(cols_names, colnames(data))
  cols_idx <- cols_idx[! is.na(cols_idx)]
  if (length(cols_idx)) {
    cols_idx <- cols_idx[unlist(lapply(cols_idx, function(x) {
      ! all(is.na(data[, x]))
    }))]
    data_names <- data[, cols_idx, drop = FALSE]
    sapply(seq_len(nrow(data_names)), function(nr) {
      entries <- lapply(data_names[nr, ! is.na(data_names[nr, ])], as.character)
      do.call(paste, as.list(c(entries, sep = "-")))
    })
  } else {
    warning("no recognized columns for entry id")
    paste0("entry", seq_len(nrow(data)))
  }
}

#' Create Row Names for STR Data
#'
#' Create unique rownames for the given data frame, using whichever STR-related
#' metadata columns are available.
#'
#' @param data STR data frame, such as produced by [prepare_dataset] or a
#'   summary from [summarize_dataset].
#'
#' @seealso [order_entries]
#'
#' @return vector of unique row names
#' @md
make_rownames <- function(data) {
  make.unique(make_entry_id(data))
}

#' Define Ordering for STR Data
#'
#' Create a standardized ordering vector (as per [base::order]) for the
#' rows of the given data frame, using whichever STR-related metadata columns
#' are available.  Sample identifiers are treated as integers primarily but then
#' resolved further by character sorting.  Note that the identification-related
#' report code relies on this to order lowest-distance genotype matches first.
#'
#' @seealso [make_rownames]
#'
#' @param data STR data frame, such as produced by [prepare_dataset] or a
#'   summary from [summarize_dataset].
#'
#' @return integer vector of new row ordering
#' @md
order_entries <- function(data) {
  items <- list(data$Locus,
                as.integer(gsub("[^0-9]+", "", data$Sample)),
                data$Sample,
                data$Dataset,
                data$Replicate,
                data$Distance,
                data$Name)
  items <- items[sapply(items, length) != 0]
  do.call(order, items)
}


# Alleles -----------------------------------------------------------------


#' Create Short Allele Names
#'
#' Autogenerate short names for sequences using sequence length and content.
#'
#' @param data character vector of sequences
#' @param hash_len number of characters of alphanumeric hash to include as a
#'   suffix.  (See [openssl::md5].)
#'
#' @return vector of short names
#' @md
make_allele_name <- function(data, hash_len = cfg("allele_suffix_len")) {
  txt <- if (is.character(data)) {
    if (hash_len > 0) {
      paste(nchar(data),
            substr(openssl::md5(data), 1, hash_len),
            sep = "-")
    } else {
      nchar(data)
    }
  } else {
    as.character(data)
  }
  txt[is.na(data)] <- NA
  txt
}

# Defines a standard order for allele names.
order_alleles <- function(nms) {
  ints <- as.integer(gsub("[^0-9]+.*", "", nms))
  order(ints, nms)
}

#' Name allele sequences in genotype data frame
#'
#' Add Allele1Name and Allele2Name columns matching Allele1Seq and Allele2Seq in
#' the given data frame.  Names from the given known_alleles data frame will be
#' used for recognized sequences.
#'
#' @param data data frame containing Allele1Seq and Allele2Seq columns such as
#'   the first list item produced by [analyze_dataset].
#' @param known_alleles data frame of custom allele names as defined for
#'   [load_allele_names].  if NULL only automatically generated names will be
#'   used.
#' @param name_args list of additional arguments to [make_allele_name].
#'
#' @return data frame provided with Allele1Name and Allele2Name columns added
#' @md
name_alleles_in_table <- function(
    data, known_alleles = NULL, ...) {
  # Make names for given seqs, using existing names where available.
  nm <- function(seqs) {
    nms <- do.call(make_allele_name, c(list(data = seqs), alist(...)))
    if (! is.null(known_alleles)) {
      idx <- match(seqs, known_alleles$Seq)
      nms <- ifelse(is.na(idx),
                    nms,
                    as.character(known_alleles$Name[idx]))
    }
    nms
  }
  # Name all of the called alleles across entries
  data$Allele1Name <- nm(data$Allele1Seq)
  data$Allele2Name <- nm(data$Allele2Seq)
  data
}


# Other -------------------------------------------------------------------

#' Remove shared path from file paths
#'
#' For the given character vector of file paths, create a modified version with
#' any common prefix path removed.  Forward slashes are used as the path
#' separator on all platforms.
#'
#' @param fps_full character vector of file paths.
#'
#' @return character vector of same length as input, with any common directory
#'   structure trimmed off.
#' @md
remove_shared_root_dir <- function(fps_full) {
  fps <- gsub("\\\\", "/", fps_full)
  fps <- normalizePath(fps, mustWork = FALSE, winslash = "/")
  chunks <- lapply(strsplit(fps, "/"), function(segs) segs[segs != ""])
  minlen <- min(sapply(chunks, length))
  dirs <- do.call(rbind, lapply(chunks, "[", 1:minlen))
  dirs <- apply(dirs, 2, function(d) length(unique(d)))
  diridx <- which(dirs > 1)[1]
  fps_rel <- sapply(chunks, function(chunk) {
    do.call(file.path, as.list(chunk[diridx:length(chunk)]))
  })
  names(fps_rel) <- fps_full
  fps_rel
}

# Equivalent of /dev/null for the build platform.
fp_devnull <- c(unix = "/dev/null", windows = "nul")[.Platform$OS.type] # nolint

#' Write Log Message
#'
#' Print a log message to the standard error stream.
#'
#' @param msg text to print.
#' @param col2 extra text to show at right margin; defaults to current time.
#' @param end ending to concatenate to message; defaults to newline character.
#' @md
logmsg <- function(msg, col2 = as.character(Sys.time()), end = "\n") {
  if (! cfg("verbose")) {
    return()
  }
  if (!is.null(col2)) {
    # right-justify col2, aim to fit total msg within 80 characters
    pad <- max(1, 80 - nchar(msg) - nchar(col2))
    msg <- paste0(msg, paste(rep(" ", pad),  collapse = ""), col2)
  }
  # stderr: file descriptor 2
  cat(paste0(msg, end), file = 2)
}

#' Reverse complement sequences
#'
#' Each entry in the input character vector is reversed and nucleotide
#' characters replaced with their complements (leaving any other text characters
#' unchanged).
#'
#' @param txt character vector of sequences
#' @returns character vector of reverse complements
#' @export
#' @md
revcmp <- function(txt) {
  bases <- CMP
  bases_lower <- tolower(CMP)
  names(bases_lower) <- tolower(names(bases_lower))
  bases <- c(bases, bases_lower)
  nas <- is.na(txt)
  txt[nas] <- ""
  out <- vapply(strsplit(txt, ""), function(vec) {
    vec <- rev(vec)
    out <- bases[vec]
    out[is.na(out)] <- vec[is.na(out)]
    paste(out, collapse = "")
  }, character(1))
  out[nas] <- NA
  out
}

# Does the given vector look "blank"?
# NULL: TRUE (caught by length of 0)
# NA: TRUE
# "": TRUE
# any combo of NA and "": TRUE
# all else: FALSE
is_blank <- function(vec) {
  length(vec) == 0 || all(vec %in% c(NA, ""))
}
