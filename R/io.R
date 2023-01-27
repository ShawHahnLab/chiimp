# Save and load data.

# Load Inputs -------------------------------------------------------------

# Expected column names for locus_attrs
locus_attrs_cols <- c("Locus", "LengthMin", "LengthMax", "LengthBuffer",
                      "Motif", "Primer", "ReversePrimer")

# Exepcted column names for genotypes_known
genotypes_cols <- c("Name", "Locus", "Allele1Seq", "Allele2Seq")

# Expected column names for allele_names
allele_names_cols <- c("Locus", "Name", "Seq")

# Expected column names for dataset
dataset_cols <- c("Filename", "Replicate", "Sample", "Locus")

#' Load configuration file
#'
#' Load a CSV or YAML-formatted text file of configuration options for
#' microsatellite analysis.  The [main] function loads configuration options
#' with this and sets them package-wide with [apply_config].
#'
#' Whatever entries are present in the file will be returned, but names that
#' don't match known names (see [CFG_DEFAULTS]) will be reported in a warning.
#'
#' @param fp path to configuration file.
#'
#' @return data frame of configuration options
#'
#' @examples
#' filename <- system.file("example_config.yml", package = "chiimp")
#' config <- load_config(filename)
#' # And then: full_analysis(config)
#' @export
#' @md
load_config <- function(fp) {
  yamls <- c("\\.yml$", "\\.yaml$")
  if (any(sapply(yamls, function(pat) grepl(pat, fp, ignore.case = TRUE)))) {
    load_config_yaml(fp)
  } else {
    load_config_csv(fp)
  }
}

# spec for any config table
CONFIG_TABLE_SPEC <- data.frame(
  Column = c(
    "Key", "Value", "Description", "Example", "Parser", "OldName", "Comments"),
  Class = "character",
  ID = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  Required = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)

check_spec <- function(obj, spec) {
  cols_req <- subset(spec, Required)$Column
  col_missing <- is.na(match(cols_req, colnames(obj)))
  if (any(col_missing)) {
    warning(paste("Missing columns in table:",
                  cols_req[col_missing]))
  }
  cols_id <- subset(spec, ID)$Column
  entries <- do.call(paste, obj[, cols_id, drop = FALSE])
  entries_repeated <- unique(entries[-match(entries, entries)])
  if (length(entries_repeated) > 0) {
    warning(paste("Repeated entries in table:", entries_repeated))
  }
  obj
}

#' Load configuration file from CSV
#'
#' @param fp path to configuration file.
#'
#' @return data frame of configuration information
#'
#' @md
load_config_csv <- function(fp) {
  cfg_table <- load_csv(fp, row.names = NULL, colClasses = "character")
  cfg_table <- check_spec(cfg_table, CONFIG_TABLE_SPEC)
  if (exists("CFG_DEFAULTS")) {
    # If we have the defaults loaded (if not, this should *be* loading the
    # defaults) check this file against those values
    config_check_keys(cfg_table)
    config_check_version(cfg_table)
  }
  cfg_table
}

#' Load configuration file from YAML
#'
#' @param fp path to configuration file.
#'
#' @return data frame of configuration information
#'
#' @md
load_config_yaml <- function(fp) {
  if (is.na(fp)) {
    config <- list()
  } else {
    text <- readChar(fp, file.info(fp)$size)
    config <- yaml::yaml.load(text)
  }
  # Flatten nested config list into table
  load_config_flatten_lists <- function(config_list, prefix = "") {
    output <- list()
    for (idx in seq_along(config_list)) {
      if (is.list(config_list[[idx]])) {
        output <- c(output, load_config_flatten_lists(
          config_list[[idx]],
          paste0(prefix, ":", names(config_list)[[idx]])))
      } else {
        output <- c(output, config_list[[idx]])
        names(output)[length(output)] <- paste0(
          prefix, ":", names(config_list)[[idx]])
      }
    }
    names(output) <- sub("^:", "", names(output))
    output
  }
  config_flat <- load_config_flatten_lists(config)
  cfg_table <- data.frame(
    OldName = names(config_flat),
    Value = as.character(config_flat),
    stringsAsFactors = FALSE)
  for (idx in seq_along(cfg_table$OldName)) {
    cfg_table$Key[idx] <- CFG_DEFAULTS$Key[
      grep(cfg_table$OldName[idx], CFG_DEFAULTS$OldName, fixed = TRUE)[1]]
  }
  cfg_table$Key[is.na(cfg_table$Key)] <- cfg_table$OldName[is.na(cfg_table$Key)]
  config_check_keys(cfg_table)
  cfg_table
}

#' Load and save tables from CSV
#'
#' Load/save a comma-separated table from/to a CSV file.  (These are generic
#' wrapper functions used by more specific loaders like [load_locus_attrs].)
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to [utils::read.table] or
#'   [utils::write.table].
#'
#' @return data frame
#'
#' @describeIn load_csv Load CSV
#'
#' @export
#' @md
load_csv <- function(fp, ...) {
  data <- utils::read.table(
    fp, header = TRUE, sep = ",", stringsAsFactors = FALSE, ...)
  if (! "row.names" %in% names(list(...))) {
    rownames(data) <- make_rownames(data)
  }
  data
}

#' @describeIn load_csv Save CSV
#' @param data data frame to save to CSV file
#' @export
save_csv <- function(data, fp, ...) {
  mkparents(fp)
  utils::write.table(
    data, file = fp, sep = ",", na = "", row.names = FALSE, ...)
  data
}

#' Load table of locus attributes
#'
#' Load a comma-separated table of locus attributes to use for analysis.  This
#' is called automatically during [full_analysis], with the data frame then used
#' by [analyze_seqs] within [analyze_dataset].
#'
#' @details
#' Columns Required:
#'   * `Locus`: Unique identifier for a given locus
#'   * `LengthMin`: Minimum known allele sequence length for this locus
#'   * `LengthMax`: Minimum known allele sequence length for this locus
#'   * `LengthBuffer`: Additional sequence length below `LengthMin` and above
#'     `LengthMax` to accept for a candidate allele
#'   * `Primer`: The forward PCR primer sequence for a given locus, used when
#'     matching sequences to loci
#'   * `ReversePrimer`: The reverse PCR primer sequence
#' @md
#'
#' @param fp_locus_attrs path to text file.
#' @param ... additional arguments passed to [load_csv].
#'
#' @return data frame of locus attributes
#'
#' @export
#'
#' @examples
#' filename <- system.file("example_locus_attrs.csv", package = "chiimp")
#' locus_attrs <- load_locus_attrs(filename)
#' @md
load_locus_attrs <- function(fp_locus_attrs, ...) {
  data <- load_csv(fp_locus_attrs, ...)
  col_missing <- is.na(match(locus_attrs_cols, colnames(data)))
  if (any(col_missing)) {
    warning(paste("Missing columns in locus_attrs table:",
                  locus_attrs_cols[col_missing]))
  }
  rows_repeated <- data$Locus[-match(data$Locus, data$Locus)]
  if (length(rows_repeated) > 0) {
    warning(paste("Repeated loci in attrs table:", rows_repeated))
  }
  rownames(data) <- make_rownames(data)
  data
}

#' Load table of known allele names
#'
#' Load a comma-separated table of allele names for use in reporting.
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to [load_csv].
#'
#' @return data frame of allele names
#'
#' @export
#' @md
load_allele_names <- function(fp, ...) {
  data <- load_csv(fp, ...)
  col_missing <- is.na(match(allele_names_cols, colnames(data)))
  if (any(col_missing)) {
    warning(paste("Missing columns in allele.names table:",
                  allele_names_cols[col_missing]))
  }
  data
}

#' Load table of genotypes
#'
#' Load a comma-separated table of genotypes, one pair of alleles per row.  This
#' information is used to compare samples to genotypes of known individuals in
#' [summarize_dataset].
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to [load_csv].
#'
#' @return data frame of genotypes
#'
#' @export
load_genotypes <- function(fp, ...) {
  data <- load_csv(fp, colClasses = "character", na.strings = "", ...)
  col_missing <- is.na(match(genotypes_cols, colnames(data)))
  if (any(col_missing)) {
    warning(paste("Missing columns in genotypes table:",
                  genotypes_cols[col_missing]))
  }
  data <- data[order_entries(data), ]
  rownames(data) <- make_rownames(data)
  data
}

#' Load table of sample attributes
#'
#' Load a comma-separated table of sample attributes for the dataset to be
#' analyzed.  Alternatively, use [prepare_dataset] to automatically
#' read sample attributes from filenames.  If more than one locus is to be
#' analyzed from a single sequencer sample (i.e., multiplexed samples), either
#' the `locusmap` argument to `prepare_dataset` can be used, or `load_dataset`
#' with an explicit mapping of loci to files.
#'
#' @details
#' Columns Required:
#' * Filename: path to each data file
#' * Replicate: identifier for repeated samples; use blanks if not applicable
#' * Sample: identifier for a given biological sample
#' * Locus: locus identifier matching that used in the locus attributes table
#'   (see [load_locus_attrs])
#' @md
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to [load_csv].
#'
#' @return data frame of sample attributes for the dataset.
#'
#' @export
#' @md
load_dataset <- function(fp, ...) {
  data <- load_csv(fp, colClasses = "character", na.strings = "", ...)
  col_missing <- is.na(match(dataset_cols, colnames(data)))
  files_missing <- ! file.exists(data$Filename)
  if (any(files_missing)) {
    logmsg(paste("WARNING: Missing", sum(files_missing), "of",
                 length(files_missing), "data files"))
  }
  if (any(col_missing)) {
    warning(paste("Missing columns in dataset table:",
                  dataset_cols[col_missing]))
  }
  # check for duplicated Sample+Replicate+Locus entries
  cts <- table(with(data, paste(Sample, Replicate, Locus, sep = "/")))
  dups <- names(cts[cts > 1])
  if (length(dups) > 0) {
    warning(paste(
      "Duplicated sample/replicate/locus entries for",
      paste(dups, collapse = ", ")))
  }
  rownames(data) <- make_rownames(data)
  data
}

#' Save table of sample attributes
#'
#' Save a comma-separated table of sample attributes.  (This is a convenience
#' function not used automatically in the analysis.)
#'
#' @param data data frame of sample attributes as produced by [prepare_dataset]
#'   or [load_dataset].
#' @param fp path to text file.
#' @param ... additional arguments passed to [save_csv].
#'
#' @export
#' @md
save_dataset <- function(data, fp, ...) {
  save_csv(data, fp, ...)
}

#' Extract Sample Attributes from Filenames
#'
#' Find files matching a pattern in a given directory, and build a data frame of
#' standard sample attributes from fields in the filenames.  Nested directory
#' structures are supported.  Alternatively, use [load_dataset] to load a
#' spreadsheet of sample attributes explicitly.  `load_dataset` can be used for
#' cases where more than one locus is to be analyzed from a single sequencer
#' sample (i.e., multiplexed samples), though the `locusmap` argument here can
#' allow automatic matching of locus names for multiplexed samples.  If the
#' directory path given does not exist or if no matching files are found, an
#' error is thrown.
#'
#' @param dp directory path to search for matching data files.
#' @param pattern regular expression to use for parsing filenames.  There should
#'   be exactly three groups in the pattern, for Replicate, Sample, and Locus.
#' @param ord integer vector giving order of the fields Replicate, Sample, and
#'   Locus in filenames.  For example, if Locus is the first field followed by
#'   Replicate and Sample, set `ord=c(3, 1, 2)`.
#' @param autorep logical allowing for automatic handling of any duplicates
#'   found, labeling them as replicates.  FALSE by default.
#' @param locusmap list of character vectors, each list item name being the
#'   locus text given in the filenames, and each vector being a set of separate
#'   locus names.  Each entry with a locus name text matching one of these list
#'   items will be replaced in the final output with several separate entries,
#'   one for each locus name in the corresponding vector.  (For example,
#'   `locusmap=list(ABCD=c("A", "B", "C", "D"))` would take a filename with
#'   "ABCD" in the locus field and split it out into four entries for the four
#'   loci.)
#'
#' @return data frame of metadata for all files found
#'
#' @export
#' @md
prepare_dataset <- function(
    dp = cfg("prep_dataset_path"),
    pattern = cfg("prep_dataset_pattern"),
    ord = cfg("prep_dataset_order"),
    autorep = cfg("prep_dataset_autorep"),
    locusmap = NULL) {
  if (! dir.exists(dp)) {
    stop(paste("ERROR: directory path for data files does not exist:",
               dp))
  }
  # get all matching filenames and extract substrings
  seq_files <- list.files(
    path = dp,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE)
  if (! length(seq_files)) {
    stop(paste("ERROR: no data files found:",
               dp))
  }
  seq_file_attrs <- stringr::str_match_all(seq_files, pattern)
  if (! all(sapply(seq_file_attrs, length) == length(ord) + 1)) {
    warning("Some filenames did not match the given pattern")
  }
  # transpose so the list is grouped by filename/attribute rather than
  # entry by entry.
  seq_file_attrs <- lapply(seq(length(ord) + 1), function(y) {
      sapply(seq_file_attrs, "[", y)
    })
  # filename as the first column, then whatever order the attrs are in
  n <- c("Filename", "Replicate", "Sample", "Locus")
  names(seq_file_attrs) <- n[c(1, 1 + ord)]
  data <- do.call(data.frame, c(seq_file_attrs, stringsAsFactors = FALSE))
  data$Filename <- seq_files
  data$Replicate <- as.integer(ifelse(data$Replicate == "",
                                      NA,
                                      data$Replicate))

  # If specified, map the locus text into multiple loci per sample (e.g.
  # multiplexed)
  if (! is_blank(locusmap)) {
    data <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
      col_locus <- match("Locus", colnames(data))
      cols <- as.list(data[i, -col_locus])
      locus <- data[i, col_locus]
      locus_new <- locusmap[[locus]]
      if (is.null(locus_new)) {
        locus_new <- locus
      }
      args <- c(cols, Locus = list(locus_new), stringsAsFactors = FALSE)
      do.call(data.frame, args)
    }))
  }

  # order by locus/sample/replicate
  data <- data[order_entries(data), ]

  # If specified, automatically number duplicates as replicates, if they don't
  # have a replicate number already.
  if (autorep) {
    data <- do.call(rbind, lapply(split(data, paste(data$Sample, data$Locus)),
           function(chunk) {
             chunk$Replicate <- ifelse(
               is.na(chunk$Replicate), seq_len(nrow(chunk)), chunk$Replicate)
             chunk
    }))
    data <- data[order_entries(data), ]
  }

  rownames(data) <- make_rownames(data)
  # complain if any replicate/sample/locus combo matches more than one entry
  if (max(table(paste(data$Replicate, data$Sample, data$Locus))) > 1) {
    warning("Some replicate/sample/locus combinations match multiple files")
  }
  return(data)
}

#' Load vector of sequences from FASTA/FASTQ file
#'
#' Load a vector of character sequences from the given path.  This is just a
#' wrapper around [dnar::read.fa] to choose the parser based on filename.  Only
#' the sequences are returned, not IDs or quality scores.
#'
#' @param fp path to sequence file
#'
#' @return vector of sequences
#'
#' @export
#' @md
load_seqs <- function(fp) {
  if (length(grep("q$", fp)) || length(grep("q.gz$", fp)))
    loadfunc <- dnar::read.fastq
  else
    loadfunc <- dnar::read.fa
  loadfunc(fp)$seq
}

# Output Saving -----------------------------------------------------------

#' Save dataset summary to text file
#'
#' Save the dataset summary produced by [analyze_dataset] to the specified file
#' path in CSV format.
#'
#' @param results_summary summary data frame as produced by [analyze_dataset].
#' @param fp output file path.
#'
#' @export
#' @md
save_results_summary <- function(results_summary, fp) {
  save_csv(results_summary, fp)
}

#' Save identified alleles to FASTA files
#'
#' Take the alleles identified by [analyze_dataset] in the summary data frame
#' and save each entry to a separate FASTA file.  Samples identified as
#' homozygous will have one sequence written rather than two.  Entries with no
#' identified alleles will be skipped.
#'
#' @param results_summary summary data frame as produced by [analyze_dataset].
#' @param dp output directory path to use for all files.
#'
#' @export
#' @md
save_allele_seqs <- function(results_summary, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  results_summary$.row <- rownames(results_summary)
  invisible(apply(results_summary, 1, function(r) {
    seqs <- r[c("Allele1Seq", "Allele2Seq")]
    # Ignore NA entries, and skip this sample if there are no sequences
    seqs <- seqs[!is.na(seqs)]
    if (length(seqs) == 0)
      return()
    # Create sequence names with a unique ID based on the sample and some other
    # stats for convenience
    seq.names <- sapply(seq_along(seqs), function(i) {
      seq.id <- paste(r[".row"], i, sep = "_")
      a_len <- r[paste0("Allele", i, "Length")]
      a_cts <- r[paste0("Allele", i, "Count")]
      seq.extra <- paste(paste0(a_len, "bp"),
                         paste0("allele/locus/total=",
                                paste(as.integer(a_cts),
                                      as.integer(r["CountLocus"]),
                                      as.integer(r["CountTotal"]),
                                      sep = "/")))
      paste(seq.id, seq.extra)
    })
    fp <- file.path(dp, paste0(r[".row"], ".fasta"))
    dnar::write.fa(seq.names, seqs, fp)
  }))
}

#' Save per-file processed data to text files
#'
#' Save each per-file data frame produced by [analyze_dataset] (via
#' [analyze_seqs]) to a separate file in the specified directory path, in CSV
#' format.  The directory structure will start at the first shared directory of
#' the input file paths. For example, if the inputs were `/data/run1/file.fastq`
#' and `/data/run2/file.fastq` there will be run1 and run2 directories inside
#' the given `dp` directory.
#'
#' @param results_file_data list of per-file data frames as produced by
#'   [analyze_dataset].
#' @param dp output directory path to use for all files.
#'
#' @export
#' @md
save_seqfile_data <- function(results_file_data, dp) {
  fps_rel <- remove_shared_root_dir(names(results_file_data))
  invisible(lapply(names(results_file_data), function(n) {
    fp_this <- fps_rel[n]
    fp <- file.path(dp, paste0(fp_this, ".csv"))
    save_csv(results_file_data[[n]], fp, quote = FALSE)
  }))
}

#' Save per-sample processed data to text files
#'
#' Save each per-sample data frame produced by [analyze_dataset] (via
#' [analyze_sample]) to a separate file in the specified directory path, in CSV
#' format.
#'
#' @param results_data list of per-sample data frames as produced by
#'   [analyze_dataset].
#' @param dp output directory path to use for all files.
#'
#' @export
#' @md
save_sample_data <- function(results_data, dp) {
  invisible(lapply(names(results_data), function(n) {
    fp <- file.path(dp, paste0(n, ".csv"))
    save_csv(results_data[[n]], fp, quote = FALSE)
  }))
}

#' Save alignments to FASTA files
#'
#' Take a list of alignments, one per locus, and save each to a separate FASTA
#' file in a specified directory.  If any of the per-locus alignment objects is
#' NA it will be skipped.  These are produced by [summarize_dataset] via
#' [align_alleles].
#'
#' @param alignments list of MSA alignment objects, such as created by
#'   [summarize_dataset] via [align_alleles].  The name of each alignment will
#'   be used for its filename.
#' @param dp output directory path.
#'
#' @export
#' @md
save_alignments <- function(alignments, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if (!is.null(alignments[[loc]])) {
      dna <- as.character(alignments[[loc]])
      fp <- file.path(dp, paste0(loc, ".fasta"))
      dnar::write.fa(names = names(dna), dna = dna, fileName = fp)
    }
  }))
}

#' Save alignment visualizations to image files
#'
#' Take a list of alignments, one per locus, and save a plot of each to a
#' separate image file in a specified directory.  If any of the per-locus
#' alignment objects is NA it will be skipped.  These are produced by
#' [summarize_dataset] via [align_alleles].
#'
#' @param alignments list of MSA alignment objects, such as created by
#'   [summarize_dataset] via [align_alleles].  The name of each alignment will
#'   be used for its filename.
#' @param dp output directory path.
#' @param image.func name of function to call for saving each image.
#' @param width integer width of image.
#' @param height integer height of image.
#' @param res integer resolution of image in PPI.
#'
#' @export
#' @md
save_alignment_images <- function(
    alignments, dp, image.func = "png",
    width = 1600, height = 1200, res = 150) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if (!is.null(alignments[[loc]])) {
      fp <- file.path(dp, paste(loc, image.func, sep = "."))
      img_call <- call(image.func,
                       fp,
                       width = width,
                       height = height,
                       res = res)
      eval(img_call)
      plot_alignment(alignments[[loc]],
                     main = paste("Alignment for Locus", loc))
      grDevices::dev.off()

    }
  }))
}

#' Save sequence histogram visualizations to image files
#'
#' Take a full results list and save a histogram (via [histogram])
#' of each sample to a separate image file in a specified directory.
#'
#' @param results list of results as created by [analyze_dataset].
#' @param dp output directory path.
#' @param image.func name of function to call for saving each image.
#' @param width integer width of image.
#' @param height integer height of image.
#' @param res integer resolution of image in PPI.
#'
#' @export
#' @md
save_histograms <- function(
    results, dp, image.func = "png", width = 1600, height = 1200, res = 150) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(results$samples), function(entry) {
    fp <- file.path(dp, paste(entry, image.func, sep = "."))
    img_call <- call(image.func, fp, width = width, height = height, res = res)
    fn <- results$summary[entry, "Filename"]
    seq_data <- results$files[[fn]]
    sample_data <- results$samples[[entry]]
    eval(img_call)
    histogram(seq_data = seq_data, sample_data = sample_data, main = entry)
    grDevices::dev.off()
  }))
}

#' Save sample distance matrix to text file
#'
#' Save the inter-sample distance matrix produced by
#' [summarize_dataset] to the specified file path in CSV format.
#'
#' @param dist_mat matrix produced by [summarize_dataset] via `make_dist_mat`.
#' @param fp output file path.
#'
#' @export
#' @md
save_dist_mat <- function(dist_mat, fp) {
  save_csv(dist_mat, fp)
}


# Helpers -----------------------------------------------------------------


# Create all parent directories as needed
mkparents <- function(fp) {
  if (!dir.exists(dirname(fp)))
    dir.create(dirname(fp), recursive = TRUE)
}
