# Save and load data.

# Expected column names for locus_attrs
locus_attrs_cols <- c("LengthMin", "LengthMax", "LengthBuffer", "Motif",
                      "Primer", "ReversePrimer")

#' Load table of locus attributes
#'
#' Load a tab-separated table of locus attributes to use for analysis.
#'
#' @param fp.locus_attrs path to text file
#'
#' @return data frame of locus attributes
#'
#' @export
load.locus_attrs <- function(fp.locus_attrs, ...) {
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
#' @param ord integer vector giving order of fields in filenames.
#'
#' @return data frame of metadata for all files found
#'
#' @export
prepare.dataset <- function(directory, pattern, ord = c(1, 2, 3)) {
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
  data$Replicate <- ifelse(data$Replicate == "", NA, data$Replicate)
  data$Locus <- ifelse(data$Locus == "", NA, data$Locus)
  rownames(data) <- make.rownames(data)
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
load.seqs <- function(fp.seqs) {
  if (length(grep("q$", fp.seqs)) || length(grep("q.gz$", fp.seqs)))
    loadfunc <- dnar::read.fastq
  else
    loadfunc <- dnar::read.fa
  loadfunc(fp.seqs)$seq
}

# create unique rownames for the given data frame, using whichever sample
# metadata columns are available.
make.rownames <- function(data) {
  cols.names <- c("Sample", "Replicate", "Locus")
  cols.idx <- match(cols.names, colnames(data))
  cols.idx <- cols.idx[!is.na(cols.idx)]
  cols.idx <- cols.idx[unlist(lapply(cols.idx, function(x) {
    !all(is.na(data[,x]))
  }) )]
  make.unique(do.call(paste, c(data[, cols.idx, drop=F], sep='-')))
}
