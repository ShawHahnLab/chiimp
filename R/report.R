# Handle report generation and plotting

# sequence histogram
# between-sample distance matrix plot

#' Render histogram of STR sample
#'
#' Given a processed STR sample, plot a histogram of counts per sequence length.
#'
#' @param samp data frame of dereplicated sequences.
#' @param stacked should bars be stacked together to add up to the total count
#'   for each length?  Otherwise they will be overlaid, so each bar height only
#'   corresponds to the count for that exact sequence.
#' @param main title of the plot.
#'
#' @export
#' @importFrom magrittr "%>%"
histogram <- function(samp, stacked=TRUE, main=NULL) {
  col.unlabeled <- "#CCCCCC" # The non-locus sequences
  col.labeled   <- "#888888" # Locus-labeled sequences
  col.flag      <- "#FF0000" # Flagged rows (presumably candidate alleles)
  lwd <- 5

  # TODO:
  #   support for showing stutter
  #   locus-specific options?
  #   un-stacked version
  #      show threshold


  if (stacked) {
    if (missing(main))
      main <- "Histogram of Sequence Lengths"
    heights <- with(samp, {
        samp %>%
        dplyr::group_by(Length) %>%
        dplyr::summarize(TotalCount = sum(Count))
      })
    # Total ranges for stacked counts
    xlim <- range(heights$Length)
    ylim <- range(heights$TotalCount)
    heights.loci <- with(samp, {
      samp %>%
      subset(MatchingLocus != "") %>%
      dplyr::group_by(Length) %>%
      dplyr::summarize(TotalCount = sum(Count))
    })

    # Plot histogram of total counts
    plot(heights$Length,
         heights$TotalCount,
         type = "h",
         xlim = xlim,
         ylim = ylim,
         main = main,
         col = col.unlabeled,
         xlab = "Length (nt)",
         ylab = "Sequence Count",
         lwd = lwd,
         lend = 1)

    # Stack counts for individual sequences for the locus-labeled cases.  If any
    # rows have the attribute "Flag" set to true, color them separately.
    for (len in unique(samp$Length)) {
      idx <- with(samp, Length == len & MatchingLocus != "")
      chunk <- samp[idx, ]
      if (nrow(chunk) < 1)
        next
      chunk <- chunk[order(chunk$Count), ]
      chunk$Cumsum <- cumsum(chunk$Count)
      for (n in 1:nrow(chunk)) {
        x <- c(len - lwd / 10,
               len - lwd / 10,
               len + lwd / 10,
               len + lwd / 10)
        y <- c(chunk[n, "Cumsum"] - chunk[n, "Count"],
               chunk[n, "Cumsum"],
               chunk[n, "Cumsum"],
               chunk[n, "Cumsum"] - chunk[n, "Count"])
        col <- col.labeled
        col.border <- "#000000"
        if (!is.null(chunk[n, "Flag"]) && chunk[n, "Flag"]) {
          col <- col.flag
          col.border <- col.flag
        }
        polygon(x, y, col = col, border = col.border)
      }
    }
  }
}
