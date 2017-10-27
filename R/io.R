# Save and load data.

# Load Inputs -------------------------------------------------------------

# Expected column names for locus_attrs
locus_attrs_cols <- c("LengthMin", "LengthMax", "LengthBuffer", "Motif",
                      "Primer", "ReversePrimer")


#' Load configuration file
#'
#' Load a YAML-formatted text file of configuration options for microsatellite
#' analysis.  This is currently just a wrapper around \code{yaml.load()}.
#'
#' @param fp.config path to configuration file.
#'
#' @return list of configuration options
#'
#' @export
load_config <- function(fp.config) {
  if (is.na(fp.config))
    return(list())
  text <- readChar(fp.config, file.info(fp.config)$size)
  yaml::yaml.load(text)
}

#' Load table of locus attributes
#'
#' Load a tab-separated table of locus attributes to use for analysis.
#'
#' @param fp.locus_attrs path to text file.
#' @param ... additional arguments passed to \code{read.table}.
#'
#' @return data frame of locus attributes
#'
#' @export
load_locus_attrs <- function(fp.locus_attrs, ...) {
  data <- read.table(fp.locus_attrs, header = T, row.names = 1, sep='\t', ...)
  col.missing <- is.na(match(locus_attrs_cols, colnames(data)))
  if (any(col.missing)) {
    warning(paste("Missing columns in locus_attrs table:",
                  locus_attrs_cols[col.missing]))
  }
  data
}

#' Load table of locus attributes
#'
#' Load a tab-separated table of locus attributes to use for analysis.
#' TODO rewrite this mess!
#'
#' @param directory location to search for matching data files.
#' @param pattern regular expression to use for parsing filenames.
#' @param ord integer vector giving order of fields (replicate/sample/locus) in
#'   filenames.
#'
#' @return data frame of metadata for all files found
#'
#' @export
prepare_dataset <- function(directory, pattern, ord = c(1, 2, 3)) {
  # get all matching filenames and extract substrings
  seq_files <- list.files(path = directory,
                          pattern = pattern,
                          full.names = TRUE,
                          recursive = TRUE,
                          include.dirs = FALSE)
  seq_file_attrs <- stringr::str_match_all(seq_files, pattern)
  # transpose so the list is grouped by filename/attribute rather than
  # entry-by-entry
  seq_file_attrs <- lapply(seq(length(seq_file_attrs[[1]])), function(y) {
      unlist(lapply(seq_file_attrs, "[", y))
    })
  # filename as the first column, then whatever order the attrs are in
  n <- c("Filename", "Replicate", "Sample", "Locus")
  names(seq_file_attrs) <- n[c(1, 1 + ord)]
  # build a data frame with filename as row name
  data <- do.call(data.frame, seq_file_attrs)
  data$Filename <- seq_files
  # order by locus/sample/replicate
  data <- data[with(data, order(Locus, Sample, Replicate)), ]
  data$Replicate <- ifelse(data$Replicate == "",
                           NA,
                           as.integer(as.character(data$Replicate)))
  #data$Locus <- ifelse(data$Locus == "", NA, data$Locus)
  rownames(data) <- make_rownames(data)
  return(data)
}

#' Load vector of sequences from FASTA/FASTQ file
#'
#' Load a vector of character sequences from the given path.  This is just a
#' wrapper around dnar to choose the parser based on filename.  Only the
#' sequences are returned, not IDs or quality scores.
#'
#' @param fp.seqs path to file
#'
#' @return vector of sequences
#'
#' @export
load_seqs <- function(fp.seqs) {
  if (length(grep("q$", fp.seqs)) || length(grep("q.gz$", fp.seqs)))
    loadfunc <- dnar::read.fastq
  else
    loadfunc <- dnar::read.fa
  loadfunc(fp.seqs)$seq
}

# Output Saving -----------------------------------------------------------

#' Save dataset summary to text file
#'
#' Save the dataset summary produced by \code{analyze_dataset} to the specified
#' file path in CSV format.
#'
#' @param results_summary summary data frame as produced by
#'   \code{analyze_dataset}.
#' @param fp output file path.
#'
#' @export
save_results_summary <- function(results_summary, fp) {
  if (!dir.exists(dirname(fp)))
    dir.create(dirname(fp), recursive = TRUE)
  write.csv(results_summary, fp, na = "")
}

#' Save identified alleles to FASTA files
#'
#' Take the alleles identified by \code{analyze_dataset} in the summary data
#' frame and save each entry to a separate FASTA file.  Samples identified as
#' homozygous will have one sequence written rather than two.  Entries with no
#' identified alleles will be skipped.
#'
#' @param results_summary summary data frame as produced by
#'   \code{analyze_dataset}.
#' @param dp output directory path to use for all files.
#'
#' @export
save_allele_seqs <- function(results_summary, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  results_summary$.row <- rownames(results_summary)
  invisible(apply(results_summary, 1, function(r) {
    seqs <- r[c("Allele1Seq", "Allele2Seq")]
    # Ignore NA entries, and skip this sample if there are no sequences
    seqs <- seqs[!is.na(seqs)]
    if (length(seqs)==0)
      return()
    # Create sequence names with a unique ID based on the sample and some other
    # stats for convenience
    seq.names <- sapply(seq_along(seqs), function(i) {
      seq.id <- paste(r[".row"], i, sep = "_")
      a.len <- r[paste0("Allele", i, "Length")]
      a.cts <- r[paste0("Allele", i, "Count")]
      seq.extra <- paste(paste0(a.len, 'bp'),
                         paste0("allele/locus/total=",
                                paste(as.integer(a.cts),
                                      as.integer(r["CountLocus"]),
                                      as.integer(r["CountTotal"]),
                                      sep = "/")))
      paste(seq.id, seq.extra)
    })
    fp <- file.path(dp, paste0(r[".row"], ".fasta"))
    dnar::write.fa(seq.names, seqs, fp)
  }))
}

#' Save per-sample processed data to text files
#'
#' Save each per-sample data frame produced by \code{analyze_dataset} to a
#' separate file in the specified directory path, in CSV format.
#'
#' @param results_data list of per-sample data frames as produced by
#'   \code{analyze_dataset}.
#' @param dp output directory path to use for all files.
#'
#' @export
save_sample_data <- function(results_data, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(results_data), function(n) {
    fp <- file.path(dp, paste0(n, '.csv'))
    write.csv(results_data[[n]], fp, na = "", quote = FALSE)
  }))
}

#' Save alignments to FASTA files
#'
#' Take a list of alignments, one per locus, and save each to a separate fasta
#' file in a specified directory.  If any of the per-locus alignment objects is
#' NA it will be skipped.
#'
#' @param alignments list of MSA alignment objects.  The name of each alignment
#'   will be used for its filename.
#' @param dp output directory path.
#'
#' @export
save_alignments <- function(alignments, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if(!is.null(alignments[[loc]])) {
      dna <- as.character(alignments[[loc]])
      fp <-file.path(dp, paste0(loc, '.fasta'))
      dnar::write.fa(names = names(dna),
                     dna = dna,
                     fileName = fp)
    }
  }))
}

#' Save alignment visualizations to image files
#'
#' Take a list of alignments, one per locus, and save a plot of each to a
#' separate image file in a specified directory.  If any of the per-locus
#' alignment objects is NA it will be skipped.
#'
#' @param alignments list of MSA alignment objects.  The name of each alignment
#'   will be used for its filename.
#' @param dp output directory path.
#' @param image.func name of function to call for saving each image.
#' @param width integer width of image.
#' @param height integer height of image.
#' @param res integer resolution of image in PPI.
#' @param ... additional arguments to \code{image.func}.
#'
#' @export
save_alignment_images <- function(alignments, dp, image.func="png",
                                  width=1600, height=1200, res=150) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if(!is.null(alignments[[loc]])) {
      fp <- file.path(dp, paste(loc, image.func, sep='.'))
      img.call <- call(image.func,
                       fp,
                       width = width,
                       height = height,
                       res = res)
      eval(img.call)
      plot_alignment(alignments[[loc]],
                     main = paste('Alignment for Locus', loc))
      dev.off()

    }
  }))
}

#' Save sample distance matrix to text file
#'
#' Save the inter-sample distance matrix produced by \code{summarize_dataset} to
#' the specified file path in CSV format.
#'
#' @param dist_mat matrix produced by \code{summarize_dataset}.
#' @param fp output file path.
#'
#' @export
save_dist_mat <- function(dist_mat, fp) {
  if (!dir.exists(dirname(fp)))
    dir.create(dirname(fp), recursive = TRUE)
  write.csv(dist_mat, fp)
}


# Misc --------------------------------------------------------------------

# create unique rownames for the given data frame, using whichever sample
# metadata columns are available.
make_rownames <- function(data) {
  cols.names <- c("Sample", "Replicate", "Locus")
  cols.idx <- match(cols.names, colnames(data))
  cols.idx <- cols.idx[!is.na(cols.idx)]
  cols.idx <- cols.idx[unlist(lapply(cols.idx, function(x) {
    !all(is.na(data[,x]))
  }) )]
  data.names <- data[, cols.idx, drop=F]
  make.unique(sapply(1:nrow(data.names), function(nr) {
    entries <- lapply((data.names[nr, !is.na(data.names[nr, ])]), as.character)
    do.call(paste, as.list(c(entries, sep='-')))
  }))
}

