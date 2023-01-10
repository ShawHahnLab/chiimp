#' Plot histogram of STR sequences
#'
#' Given processed STR sequences and optionally the filtered version, plot a
#' histogram of counts per sequence length.
#'
#' @param seq_data data frame of dereplicated sequences as created by
#'   \code{\link{analyze_seqs}}.
#' @param sample_data data frame of filtered and categorized sequences as
#'   created by \code{\link{analyze_sample}}.
#' @param main title of the plot.
#' @param xlim numeric range for x-axis.
#' @param cutoff_fraction numeric threshold for the fraction of locus-matching
#'   counts needed to call an allele.  Used to draw a horizontal line if
#'   \code{sample_data} is given.
#' @param render Should the plot be drawn to the display device?
#'
#' @return list of data frames for the sets of counts-versus-length bars drawn
#'   in the plot, split by category.
#'
#' @export
histogram <- function(seq_data,
                     sample_data = NULL,
                     main = "Number of Reads by Sequence Length",
                     xlim = range(seq_data$Length),
                     cutoff_fraction = NULL,
                     render = TRUE) {
  # TODO:
  # mark top edges of allele sequences (to handle the same-length case)
  # label bars$topknown by name? positioning may get tricky.
  bars <- str_hist_setup(seq_data, sample_data)
  if (render && nrow(seq_data) > 0) {
    if (is.null(cutoff_fraction)) {
      cutoff_fraction <- attr(sample_data, "fraction.min")
    }
    str_hist_render(bars, main, xlim, cutoff_fraction)
  }
  return(invisible(bars))
}

#' Prepare histogram data
#'
#' @param seq_data data frame of dereplicated sequences as created by
#'   \code{\link{analyze_seqs}}.
#' @param sample_data data frame of filtered and categorized sequences as
#'   created by \code{\link{analyze_sample}}.
#'
#' @return list of data frames for the sets of counts-versus-length bars drawn
#'   in the plot, split by category.
str_hist_setup <- function(seq_data, sample_data = NULL) {
  vec_to_df <- function(data, cols = c("Length", "Count")) {
    if (length(data) == 0) {
      df <- data.frame(Length = as.integer(NA), Count = as.integer(NA))[0, ]
    } else {
      df <- data.frame(as.integer(names(data)), data)
    }
    colnames(df) <- cols
    df
  }
  bars <- list()
  bars$orig <- vec_to_df(sapply(split(seq_data, seq_data$Length),
                      function(chunk) sum(chunk$Count)))

  if (!is.null(sample_data)) {
    bars$filt <- vec_to_df(sapply(split(sample_data, sample_data$Length),
                        function(chunk) sum(chunk$Count)))
    bars$topcounts <- vec_to_df(sapply(split(sample_data, sample_data$Length),
                             function(chunk) max(chunk$Count)))
    bars$topknown <- do.call(rbind, c(list(data.frame(Length = as.integer(NA),
                                               Count = as.integer(NA),
                                               SeqName = as.character(NA),
                                               stringsAsFactors = FALSE)[0, ]),
                             lapply(split(sample_data,
                                                 sample_data$Length),
                                 function(chunk) {
                                   chunk <- chunk[! is.na(chunk$SeqName), ]
                                   if (nrow(chunk) == 0)
                                     return()
                                   idx <- which(chunk$Count == max(chunk$Count))
                                   chunk[idx, c("Length", "Count", "SeqName")]
                                 })))
    bars$allele <- data.frame(Length = sample_data[sample_data$Category ==
                                                     "Allele",
                                                   "Length"],
                              Count = sample_data[sample_data$Category ==
                                                    "Allele",
                                                  "Count"])
  }
  bars <- lapply(bars, function(b) {
    rownames(b) <- NULL
    b
  })
  bars
}

#' Draw prepared histogram data
#'
#' Render prepared histogram data to the display device.
#'
#' @param bars list of data frames of counts-vs-lengths as prepared by
#' \code{\link{str_hist_setup}}.
#' @param main title of the plot.
#' @param xlim numeric range for x-axis.
#' @param cutoff_fraction numeric threshold for the fraction of locus-matching
#'   counts needed to call an allele.  Used to draw a horizontal line if
#'   \code{sample_data} is given.
str_hist_render <- function(bars, main, xlim, cutoff_fraction) {

  categories <- str_hist_setup_legend(bars)

  ylim <- range(bars$orig$Count)
  graphics::plot(
    c(),
    c(),
    main = main,
    xlab = "Length (nt)",
    ylab = "Sequence Count",
    xlim = xlim,
    ylim = ylim)

  # How wide should each bar be, in pixels, to be flush with the adjacent bars?
  lwd <- max(1, get_px_width()) # at least one pixel

  for (nm in names(bars)) {
    if (nrow(bars[[nm]]) > 0) {
      graphics::points(
        bars[[nm]]$Length,
        bars[[nm]]$Count,
        type = "h",
        col = categories[nm, "col"],
        lend = 1,
        lwd = lwd)
    }
  }

  if (! is.null(bars$filt) && nrow(bars$filt) > 0) {
    # Draw threshold
    cutoff <- (cutoff_fraction * sum(bars$filt$Count))[1]
    if (! is.na(cutoff)) {
      categories["threshold", "Render"] <- TRUE
      xlim_view <- graphics::par("usr")[1:2]
      xlim_view[1] <- floor(xlim_view[1])
      xlim_view[2] <- ceiling(xlim_view[2])
      graphics::points(xlim_view[1]:xlim_view[2],
             rep(cutoff, diff(xlim_view) + 1), pch = ".")
    }
    # Draw domain of sample data
    ymax <- max(bars$orig$Count)
    xlim_filt <- range(as.integer(bars$filt$Length))
    graphics::polygon(
      x = rep(xlim_filt, each = 2),
      y = c(0, ymax, ymax, 0),
      col = "#0000001E",
      border = NA)
  }

  # Draw legend
  leg <- categories[categories$Render, c("legend", "col", "pch", "lty")]
  filt <- sapply(leg, function(x) all(is.na(x)))
  leg <- leg[! filt]
  do.call(graphics::legend, c(list(x = "topright", bty = "n"), leg))
}

#' Setup display attributes for STR histogram
#'
#' Create data frame of plot attributes to use in STR histogram.
#'
#' @param bars data frames of counts-vs-lengths as prepared by
#' \code{\link{str_hist_setup}}.
#'
#' @return data frame of plot attributes by category.
str_hist_setup_legend <- function(bars) {
  categories <- data.frame(
    Name = c("orig", "filt", "topcounts", "topknown", "allele", "threshold"),
    Render = FALSE,
    legend = c("Original",
               "Filtered",
               "Filtered, top unique Seq.",
               "Filtered, top unique Seq., Known",
               "Called Alleles",
               "Unique Seq. Threshold"),
    col = c("black", "gray", "pink", "blue", "red", "#00000080"),
    pch = c(15, 15, 15, 15, 15, NA),
    lty = c(NA, NA, NA, NA, NA, 1),
    stringsAsFactors = FALSE)
  categories$Render[seq_len(length(bars))] <- TRUE
  rownames(categories) <- categories$Name
  categories
}

#' Get plot unit width
#'
#' How wide, in pixels, is an increment of one plot unit for the current display
#' device?
get_px_width <- function() {
  # https://stackoverflow.com/questions/17213293/how-to-get-r-plot-window-size
  px_width_fig <- grDevices::dev.size("px")[1] # width of whole figure in pixels
  px_width_plt <- diff(graphics::par("plt")[1:2] * px_width_fig) # plt region
  width_plt <- diff(graphics::par("usr")[1:2]) # plot region (in plot units)
  step_width <- px_width_plt / width_plt # pixels per plot unit increment
  step_width
}
