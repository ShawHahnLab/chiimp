# Save and load data.

#' Load table of locus attributes
#'
#' Load a tab-separated table of locus attributes to use for analysis.
#'
#' @param fp.locus_attrs path to text file
#'
#' @return data frame of locus attributes
#'
#' @export
load.locus_attrs <- function(fp.locus_attrs) {
  read.table(fp.locus_attrs, sep = "\t", header = T, row.names = 1)
}

load.seqs <- function(fp.seqs) {
  if (length(grep("q$", fp.seqs)) || length(grep("q.gz$", fp.seqs)))
    loadfunc <- dnar::read.fastq
  else
    loadfunc <- dnar::read.fa
  loadfunc(fp.seqs)$seq
}

# TODO rewrite this mess
prepare.dataset <- function(directory, pattern, ord = c(1, 2, 3)) {
  # get all matching filenames and extract substrings
  seq_files <- list.files(directory, pattern, full.names = TRUE)
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
  #row.names(data) <- seq_file_attrs[[1]]
  # order by locus/sample/replicate
  data <- data[with(data, order(Locus, Sample, Replicate)), ]
  data$Replicate <- ifelse(data$Replicate == "", NA, data$Replicate)
  rownames(data) <- NULL
  return(data)
}
