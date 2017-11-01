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

# Autogenerate a short name for a sequence using sequence length and content.
make_allele_name <- function(data, hash.len=6) {
  if (is.character(data)) {
    if(hash.len > 0) {
      paste(nchar(data),
            substr(openssl::md5(data), 1, hash.len),
            sep = "-")
    } else {
      nchar(data)
    }
  } else {
    as.character(data)
  }
}

# Return an order for the rows of the given data frame, using whichever sample
# metadata columns are available.  Treat sample identifiers as integers
# primarily but resolve further by character sorting.
order_entries <- function(data) {
  items <- list(data$Locus,
                as.integer(gsub("[^0-9]+", "", data$Sample)),
                data$Sample,
                data$Replicate)
  items <- items[!sapply(items, is.null)]
  do.call(order, items)
}
