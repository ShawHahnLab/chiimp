# Misc utility functions and variables used by the others.

#' Create Row Names for STR Data
#'
#' Create unique rownames for the given data frame, using whichever STR-related
#' metadata columns are available.
#'
#' @param data STR data frame, such as produced by \code{\link{prepare_dataset}}
#'   or \code{summary} from \code{\link{summarize_dataset}}.
#'
#' @seealso \code{\link{order_entries}}
#'
#' @return vector of unique row names
make_rownames <- function(data) {
  cols.names <- c("Dataset", "Sample", "Name", "Replicate", "Locus")
  cols.idx <- match(cols.names, colnames(data))
  cols.idx <- cols.idx[!is.na(cols.idx)]
  cols.idx <- cols.idx[unlist(lapply(cols.idx, function(x) {
    !all(is.na(data[, x]))
  }) )]
  data.names <- data[, cols.idx, drop = F]
  make.unique(sapply(1:nrow(data.names), function(nr) {
    entries <- lapply(data.names[nr, !is.na(data.names[nr, ])], as.character)
    do.call(paste, as.list(c(entries, sep = "-")))
  }))
}

#' Define Ordering for STR Data
#'
#' Create a standardized ordering vector (as per \code{\link{order}}) for the
#' rows of the given data frame, using whichever STR-related metadata columns
#' are available.  Sample identifiers are treated as integers primarily but then
#' resolved further by character sorting.
#'
#' @seealso \code{\link{make_rownames}}
#'
#' @param data STR data frame, such as produced by \code{\link{prepare_dataset}}
#'   or \code{summary} from \code{\link{summarize_dataset}}.
#'
#' @return integer vector of new row ordering
order_entries <- function(data) {
  items <- list(data$Locus,
                as.integer(gsub("[^0-9]+", "", data$Sample)),
                data$Name,
                data$Sample,
                data$Dataset,
                data$Replicate,
                data$Distance)
  items <- items[sapply(items, length) != 0]
  do.call(order, items)
}

#' Create Short Allele Names
#'
#' Autogenerate short names for sequences using sequence length and content.
#'
#' @param data character vector of sequences
#' @param hash.len number of characters of alphanumeric hash to include as a
#'   suffix.  (See \code{\link[openssl]{md5}}.)
#'
#' @return vector of short names
make_allele_name <- function(data, hash.len=6) {
  if (is.character(data)) {
    if (hash.len > 0) {
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

# Equivalent of /dev/null for the build platform.
fp_devnull <- c(unix = "/dev/null", windows = "nul")[.Platform$OS.type] # nolint
