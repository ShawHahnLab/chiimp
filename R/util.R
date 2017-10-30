# Misc utility functions used by the others.

# create unique rownames for the given data frame, using whichever sample
# metadata columns are available.
make_rownames <- function(data) {
  cols.names <- c("Sample", "Name", "Replicate", "Locus")
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
