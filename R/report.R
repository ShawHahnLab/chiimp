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
#' @param locus.name name of the locus to match alleles for.  If unspecified the
#'   locus with the highest matched counts will be used.
#' @param sample.summary summary data frame as prepared by
#'   \code{summarize.sample}
#'
#' @export
#' @importFrom magrittr "%>%"
histogram <- function(samp, stacked=TRUE, main=NULL, locus.name=NULL, sample.summary=NULL) {
  col.unlabeled   <- "#DDDDDD" # The non-locus sequences
  # Locus-labeled sequences
  col.labeled     <- "#999999"
  border.labled   <- "#000000"
  # Locus-labeled and matching all allele conditions
  col.filtered    <- "#FF6666"
  border.filtered <- "#990000"
  # Sequences matching alleles in sample.summary
  col.allele      <- "#FF0000"
  border.allele   <- "#FF0000"
  lwd <- 5

  if (missing(locus.name)) {
    cts <- subset(sample.data, MatchingLocus != "") %>%
      dplyr::group_by(MatchingLocus) %>%
      dplyr::summarize(Count=sum(Count))
    cts <- cts[order(cts$Count, decreasing = T), ]
    locus.name <- as.character(cts[[1, 'MatchingLocus']])
  }

  # TODO:
  #   support for showing stutter
  #   locus-specific options?
  #   un-stacked version
  #      show threshold

  if (stacked) {
    if (missing(main))
      main <- "Histogram of Sequence Lengths"
    # total counts at each length
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

    # Plot histogram of total counts, gray in background
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
    # rows match specific conditions, color those rectangles separately.
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
        # Is the sequence matching all criteria for a candidate allele?
        allele_match <- allele.match(chunk[n, ], locus.name)
        if(is.na(allele_match)) allele_match <- F
        # Is the sequence an exact match for an identified allele?
        is_allele <- chunk[n, 'Seq'] %in%
          sample.summary[c('Allele1Seq', 'Allele2Seq')]
        # Color each rectangle according to its category as defined by the
        # above.  The default will be the locus-labeled colors.
        col <- col.labeled
        col.border <- border.labled
        if (is_allele) {
          col <- col.allele
          col.border <- border.allele
        } else if (allele_match) {
          col <- col.filtered
          col.border <- border.filtered
        }
        polygon(x, y, col = col, border = col.border)
      }
    }
  }
}

report.genotypes <- function(dataset.summary, ...) {
  tbl <- summarize.genotypes(dataset.summary, ...)
  # If we have no replicates drop that column
  if (all(is.na(tbl$Replicate)))
    tbl <- tbl[, -2]
  # Blank out NA values
  tbl[is.na(tbl)] <- ''
  tbl
}

report.attribute <- function(dataset.summary, attrib, ...) {
  tbl <- summarize.attribute(dataset.summary, attrib, ...)
  # If we have no replicates drop that column
  if (all(is.na(tbl$Replicate)))
    tbl <- tbl[, -2]
  tbl
}

# A skewed 0 -> 1 scale for color-coding distance tables
make.dist_scale <- function(n) {
  ((0:n)/n)**(1/3)
}

plot.dist_mat <- function(dist_mat, num.alleles=max(dist_mat),
                          dist.display_thresh=round(num.alleles*2/3)) {
  labels <- matrix(character(length(dist_mat)), nrow=nrow(dist_mat))
  idx <- dist_mat <= dist.display_thresh
  labels[idx] <- dist_mat[idx]
  #idx.diag <- (0:(nrow(dist_mat)-1))*(nrow(dist_mat)+1) + 1
  #labels[idx.diag] <- ''
  diag(labels) <- ''
  dist_scale <- make.dist_scale(num.alleles)
  color <- rgb(red=1, green=dist_scale, blue=dist_scale)

  pheatmap::pheatmap(dist_mat,
           color = color,
           display_numbers = labels,
           treeheight_row = 0,
           breaks = 0:num.alleles,
           clustering_distance_rows = as.dist(dist_mat),
           clustering_distance_cols = as.dist(dist_mat))
}
