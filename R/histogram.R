# TODO:
# test with empty seq_data
# test with just seq_data
# test with seq_data + sample_data
# test with seq_data + empty sample_data

str_hist <- function(seq_data,
                     sample_data = NULL,
                     main = "Number of Reads by Sequence Length",
                     xlim = range(seq_data$Length),
                     cutoff_fraction = NULL,
                     render = TRUE) {
  bars <- str_hist_setup(seq_data, sample_data)
  if (render) {
    if (is.null(cutoff_fraction)) {
      cutoff_fraction <- attr(sample_data, "fraction.min")
    }
    str_hist_render(bars, main, xlim, cutoff_fraction)
  }
  return(invisible(bars))
}

str_hist_setup <- function(seq_data, sample_data = NULL) {
  bars <- list()
  bars$orig <- sapply(split(seq_data, seq_data$Length),
                      function(chunk) sum(chunk$Count))

  if (!is.null(sample_data)) {
    bars$filt <- sapply(split(sample_data, sample_data$Length),
                        function(chunk) sum(chunk$Count))
    bars$topcounts <- sapply(split(sample_data, sample_data$Length),
                             function(chunk) max(chunk$Count))
    bars$topknown <- unlist(sapply(split(sample_data, sample_data$Length),
                                   function(chunk) {
                                     chunk <- subset(chunk, ! is.na(SeqName))
                                     if (nrow(chunk) == 0)
                                       return()
                                     max(chunk$Count)
                                   }))
    bars$allele <- sample_data[sample_data$Category == "Allele",
                               "Count"]
    names(bars$allele) <- sample_data[sample_data$Category == "Allele",
                                      "Length"]
  }
  bars
}

str_hist_render <- function(bars, main, xlim, cutoff_fraction) {

  categories <- str_hist_setup_legend(bars)

  ylim <- range(bars$orig)
  graphics::plot(c(),
                 c(),
                 main = main,
                 xlab = "Length (nt)",
                 ylab = "Sequence Count",
                 xlim = xlim,
                 ylim = ylim)

  # How wide should each bar be, in pixels, to be flush with the adjacent bars?
  lwd <- max(1, get_px_width()) # at least one pixel

  for (nm in names(bars)) {
    graphics::points(names(bars[[nm]]),
                     bars[[nm]],
                     type = "h",
                     col = categories[nm, "col"],
                     lend = 1,
                     lwd = lwd)
  }

  if (! is.null(bars$filt)) {
    # Draw threshold
    cutoff <- (cutoff_fraction * sum(bars$filt))[1]
    if (! is.na(cutoff)) {
      categories["threshold", "Render"] <- TRUE
      points(1:xlim[2], rep(cutoff, xlim[2]), pch='.')
    }
    # Draw domain of sample data
    ymax <- max(bars$orig)
    xlim_filt <- range(as.integer(names(bars$filt)))
    graphics::polygon(x = rep(xlim_filt, each = 2),
                      y = c(0, ymax, ymax, 0),
                      col = "#0000001E",
                      border = NA)
  }

  # Draw legend
  leg <- subset(categories, Render)[, c("legend", "col", "pch", "lty")]
  filt <- sapply(leg, function(x) all(is.na(x)))
  leg <- leg[! filt]
  do.call(graphics::legend,
          c(list(x = "topright"),
            leg))
}

# default plot category attributes
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
    stringsAsFactors = FALSE
  )
  categories$Render[1:length(bars)] <- TRUE
  rownames(categories) <- categories$Name
  categories
}

# how wide, in pixels, is an increment of one plot unit?
# https://stackoverflow.com/questions/17213293/how-to-get-r-plot-window-size
get_px_width <- function() {
  px_width_fig <- dev.size("px")[1] # width of whole figure in pixels
  px_width_plt <- diff(par('plt')[1:2] * px_width_fig) # just the plot region
  width_plt <- diff(par('usr')[1:2]) # width of plot region in plot units
  step_width <- px_width_plt / width_plt # pixels per plot unit increment
  step_width
}
