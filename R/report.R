# Handle summary table generation and plotting.  No real processing should
# happen here.

# Tables ------------------------------------------------------------------

# For each pair of alleles in the given data frame, reorder in standard order,
# duplicate for homozygous, and fill in blank for remaining NAs.
normalize_alleles <- function(data) {
  as.data.frame(t(apply(data, 1, function(pair) {
    ord <- order_alleles(pair)
    pair <- pair[ord]
    if (is.na(pair[1])) {
      pair[1] <- ""
    }
    pair[is.na(pair)] <- pair[1]
    pair
  }
  )), stringsAsFactors = FALSE)
}

#' Wide table of allele names vs loci
#'
#' Allele pairs are shown in a standardized order with homozygous entries shown
#' twice.
#'
#' @param data data frame containing Allele1Name and Allele2Name colums such as
#'   the first list item produced by \code{\link{analyze_dataset}}.  If allele
#'   names are not yet present call \code{\link{name_alleles_in_table}}.
#' @param extra_cols names or index values of additional columns from input data
#'   frame to be kept in output data frame.  These should be consistent across
#'   loci for a given entry.
#'
#' @return wide format data frame with sample entries on rows and loci on
#'   columns.  An ID column will label sample entries by whichever columns were
#'   provided in the input (see \code{\link{make_entry_id}}).
#'
#' @export
tabulate_allele_names <- function(data, extra_cols=NULL) {
  # Order and replicate (for homozygous) the allele names
  nms <- normalize_alleles(data[, c("Allele1Name", "Allele2Name")])
  # Create unique (aside from Locus) identifiers for each entry
  id <- make_entry_id(data[, -match("Locus", colnames(data))])
  # Our normalized and ordered long-format data frame to be reshaped.
  long <- data.frame(ID = id,
                     Locus = data$Locus,
                     data[, extra_cols, drop = FALSE],
                     nms,
                     stringsAsFactors = FALSE)
  long <- long[order_entries(long), ]
  # Switch to wide format, putting the allele names per locus across columns
  # (along with ID and whatever extra_cols were given).
  tbl <- stats::reshape(long, v.names = c("V1", "V2"),
                        idvar = "ID",
                        timevar = "Locus",
                        direction = "wide")
  # Fix row and colunn names
  rownames(tbl) <- tbl$ID
  loci <- data$Locus
  if (is.factor(loci)) {
    loci <- levels(droplevels(data$Locus))
  } else {
    loci <- unique(loci)
  }
  allele_cols <- paste(rep(loci, each = 2),
                       c(1, 2),
                       sep = "_")
  colnames(tbl) <- c("ID", extra_cols, allele_cols)
  tbl
}

#' Create genotype summary table
#'
#' Report the genotypes present in a processed dataset in a concise data frame.
#' This will arrange the allele names into a wide-format table with unique
#' samples on rows and loci on columns, do some automatic cleanup on the
#' columns, and show closest-matching individuals per entry, if given.
#'
#' @param results list of results data as produced by \code{analyze_dataset}.
#' @param na.replicates text to replace NA entries with for the Replicates
#'     column.
#' @param closest list of closest matches as produced by
#'   \code{\link{find_closest_matches}}.
#'
#' @return data frame showing summary of genotypes.
#' @export
report_genotypes <- function(results,
                             na.replicates="",
                             closest=NULL) {
  tbl <- tabulate_allele_names(data = results$summary,
                               extra_cols = c("Sample", "Replicate"))
  tbl <- tbl[, -match("ID", colnames(tbl))]

  # If a closest-matches list was given, add columns using that-- but only if
  # there's exactly one.
  if (! is.null(closest)) {
    idents <- do.call(rbind, lapply(closest, function(who) {
      if (length(who) == 1) {
        data.frame(Distance = who, Name = names(who), stringsAsFactors = FALSE)
      } else {
        data.frame(Distance = NA, Name = NA)
      }
    }))
    tbl <- cbind(tbl, idents)
  }

  # If we have no replicates drop that column
  if (all(is.na(tbl$Replicate)))
    tbl <- tbl[, -2]
  else
    tbl$Replicate[is.na(tbl$Replicate)] <- na.replicates
  # Blank out any remaining NA values
  tbl[is.na(tbl)] <- ""

  tbl
}

#' Create identification summary table
#'
#' Report the genotypes present in a processed dataset, paired with close
#' matches to known individuals, converting sequences to short names.  This is a
#' more detailed view than given by \code{\link{report_genotypes}}.
#'
#' @param results results list as produced by \code{\link{summarize_dataset}}.
#' @param closest list of closest matches as produced by
#'   \code{\link{find_closest_matches}}.
#' @param na.replicates text to replace NA entries with for the Replicates
#'     column.
#'
#' @return data frame showing summary of sample genotypes with interleved
#'   genotypes for similar known individuals.
#'
#' @export
report_idents <- function(results,
                          closest,
                          na.replicates="") {
  # Take the known genotypes table, but keep only those entries relevant to the
  # current loci and match levels.
  gt <- subset(results$genotypes.known,
               Locus %in% levels(results$summary$Locus))
  gt$Locus <- factor(gt$Locus, levels = levels(results$summary$Locus))
  # Add Allele1Name, Allele2Name columns
  gt <- name_alleles_in_table(data = gt, known_alleles = results$allele.names)
  # Create wide-format allele name table for experimental and known values
  results_summary <- results$summary
  tbl_exp <- tabulate_allele_names(data = results_summary,
                    extra_cols = c("Sample", "Replicate"))

  tbl_known <- tabulate_allele_names(data = gt,
                    extra_cols = c("Name"))


  # Now, to group each row in tbl_exp with each set of nearby individuals in
  # tbl_known.

  # Create tbl_closest as a full genotype table of those known entries (possibly
  # multiple per sample) that are closest to each experimental sample.
  tbl_closest <- do.call(rbind, lapply(names(closest), function(nm) {
    if (length(closest[[nm]]) == 0)
      return(NULL)
    idx_known <- match(names(closest[[nm]]), tbl_known$Name)
    idx_exp <- match(nm, tbl_exp$ID)
    cbind(Sample = tbl_exp[idx_exp, "Sample"],
          Replicate = tbl_exp[idx_exp, "Replicate"],
          tbl_known[idx_known, ], # cols ID, Name, loci
          Distance = closest[[nm]])
  }))

  # Interleave with the samples themselves and group each set in a report table.
  # This will order by sample by distance, with observations first, before known
  # genotypes.
  tbl_combo <- tbl_exp
  tbl_combo$Distance <- -1
  tbl_combo$Name <- ""
  tbl_combo <- rbind(tbl_combo, tbl_closest)
  tbl_combo <- tbl_combo[order_entries(tbl_combo), ]
  tbl_combo <- within(tbl_combo, Distance[Distance < 0] <- NA)

  # Drop ID column
  tbl_combo <- tbl_combo[, -match("ID", colnames(tbl_combo))]
  # If we have no replicates drop that column
  if (all(is.na(tbl_combo$Replicate)))
    tbl_combo <- tbl_combo[, -match("Replicate", colnames(tbl_combo))]
  else
    tbl_combo$Replicate[is.na(tbl_combo$Replicate)] <- na.replicates
  # Blank out any remaining NA values
  tbl_combo[is.na(tbl_combo)] <- ""

  return(tbl_combo)
}

# Plots -------------------------------------------------------------------

#' Plot basic histogram of STR sample
#'
#' Given a processed STR sample, plot a histogram of counts per sequence length.
#' If selected, bars for the counts of sequences matching a selected locus, and
#' the counts of sequences representing called alleles, will be overlayed.
#'
#' @param samp data frame of dereplicated sequences.
#' @param main title of the plot.
#' @param locus.name name of the locus to match alleles for.  If unspecified the
#'   locus with the highest matched counts will be used.  To disable the
#'   matching entirely use \code{NA}.
#' @param sample.summary summary data frame as prepared by
#'   \code{summarize_sample}.  Used to label the called alleles.
#' @param cutoff_fraction numeric threshold for the fraction of locus-matching
#'   counts needed to call an allele.  To disable display use \code{NA}.
#' @param xlim numeric range for x-axis.
#'
#' @export
#' @importFrom magrittr "%>%"
histogram <- function(samp,
                      main="Number of Reads by Sequence Length",
                      locus.name=NULL,
                      sample.summary=NULL,
                      cutoff_fraction=0.05,
                      xlim=range(samp$Length)) {

  ## Define colors and other plot parameters
  hist_colors <- c(bar_unlabeled = "#000000FF", # non-locus sequences
                   bar_filtered  = "#888888FF", # matching allele conditions
                   bar_topcounts = "#FFAAAAFF", # , but just seq w/ max counts
                   bar_allele    = "#FF0000FF", # matching alleles in summary
                   line_cutoff   = "#00000080", # threshold for allele calls
                   rect_region   = "#0000001E") # region for filtered sequences

  lwd <- 5

  ## If there were no sequences at all, don't bother plotting.
  if (nrow(samp) == 0)
    return(NULL)

  ## Group sequences by length
  heights <- with(samp, {
    samp %>%
    dplyr::group_by(Length) %>%
    dplyr::summarize(TotalCount = sum(Count))
  })
  ymax <- max(heights$TotalCount)

  # Plot bars for all counts
  ylim <- range(heights$TotalCount)
  graphics::plot(heights$Length,
       heights$TotalCount,
       type = "h",
       xlim = xlim,
       ylim = ylim,
       main = main,
       col = hist_colors["bar_unlabeled"],
       xlab = "Length (nt)",
       ylab = "Sequence Count",
       lwd = lwd,
       lend = 1)

  # If no locus name was given, take whatever locus showed the highest counts.
  if (missing(locus.name)) {
    cts <- with(samp, {
      subset(samp, MatchingLocus != "") %>%
        dplyr::group_by(MatchingLocus) %>%
        dplyr::summarize(Count = sum(Count))
    })
    cts <- cts[order(cts$Count, decreasing = T), ]
    locus.name <- if (nrow(cts) > 1) {
      cts[[1, "MatchingLocus"]]
    } else {
      NA
    }
  }

  # Just the sequences matching all locus conditions, if there are any.
  samp.filt <- samp[full_locus_match(samp, locus.name), ]
  if (nrow(samp.filt) > 0) {
    heights.filt <- with(samp.filt, {
      samp.filt %>%
        dplyr::group_by(Length) %>%
        dplyr::summarize(TotalCount = sum(Count))
    })
    graphics::points(heights.filt$Length,
           heights.filt$TotalCount,
           type = "h",
           col = hist_colors["bar_filtered"],
           lwd = lwd,
           lend = 1)
    # Draw bars for the highest unique sequence at each length
    idx <- match(samp.filt$Length, samp.filt$Length)
    heights.filt.top <- samp.filt[idx, ]
    graphics::points(heights.filt.top$Length,
                     heights.filt.top$Count,
                     type = "h",
                     col = hist_colors["bar_topcounts"],
                     lwd = lwd,
                     lend = 1)
    # Shade the domain of the filtered data in gray
    xlim.filt <- range(samp.filt$Length)
    graphics::polygon(x = rep(xlim.filt, each = 2),
            y = c(0, ymax, ymax, 0),
            col = hist_colors["rect_region"],
            border = NA)
    # Draw a line to mark the cutoff value for what's considered a prominent
    # count
    cutoff <- cutoff_fraction * sum(samp.filt$Count)
    graphics::abline(h = cutoff, col = hist_colors["line_cutoff"])
  }

  # Draw bars for the exact allele sequences identified
  if (!is.null(sample.summary)) {
    pts.x <- c(sample.summary$Allele1Length, sample.summary$Allele2Length)
    pts.y <- c(sample.summary$Allele1Count, sample.summary$Allele2Count)
    pts.x <- pts.x[!is.na(pts.x)]
    pts.y <- pts.y[!is.na(pts.y)]
    if (length(pts.x) > 0)
      graphics::points(pts.x,
             pts.y,
             type = "h",
             col = hist_colors["bar_allele"],
             lwd = lwd,
             lend = 1)
  }

  # Legend
  graphics::legend(x = "topright", bty = "n",
         legend = c("Original",
                    "Filtered",
                    "Filtered, top unique Seq.",
                    "Called Alleles",
                    "Unique Seq. Threshold"),
         col = hist_colors[1:5],
         pch = c(15, 15, 15, 15, NA),
         lty = c(NA, NA, NA, NA, 1))
}

#' Plot advanced histogram of STR sample
#'
#' Given a processed STR sample, plot a histogram of counts per unique sequence.
#' This is a more complicated version of \code{histogram}.
#'
#' @param samp data frame of dereplicated sequences.
#' @param stacked should bars be stacked together to add up to the total count
#'   for each length?  Otherwise they will be overlaid, so each bar height only
#'   corresponds to the count for that exact sequence.
#' @param main title of the plot.
#' @param locus.name name of the locus to match alleles for.  If unspecified the
#'   locus with the highest matched counts will be used.
#' @param sample.summary summary data frame as prepared by
#'   \code{summarize_sample}
#'
#' @export
#' @importFrom magrittr "%>%"
histogram2 <- function(samp,
                      stacked=TRUE,
                      main=NULL,
                      locus.name=NULL,
                      sample.summary=NULL) {
  col.unlabeled   <- "#DDDDDD" # The non-locus sequences
  # Locus-labeled sequences
  col.labeled     <- "#999999"
  border.labled   <- "#000000"
  # Locus-labeled and matching all allele conditions
  col.filtered    <- "#FFDDDD"
  border.filtered <- "#990000"
  # Sequences matching alleles in sample.summary
  col.allele      <- "#FF0000"
  border.allele   <- "#FF0000"
  lwd <- 5

  if (missing(locus.name)) {
    cts <- with(samp, {
      subset(samp, MatchingLocus != "") %>%
        dplyr::group_by(MatchingLocus) %>%
        dplyr::summarize(Count = sum(Count))
    })
    cts <- cts[order(cts$Count, decreasing = T), ]
    locus.name <- cts[[1, "MatchingLocus"]]
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
    graphics::plot(heights$Length,
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
        am <- full_locus_match(chunk[n, ], locus.name)
        if (is.na(am)) am <- F
        # Is the sequence an exact match for an identified allele?
        is_allele <- chunk[n, "Seq"] %in%
          unlist(sample.summary[c("Allele1Seq", "Allele2Seq")])
        # Color each rectangle according to its category as defined by the
        # above.  The default will be the locus-labeled colors.
        col <- col.labeled
        col.border <- border.labled
        if (is_allele) {
          col <- col.allele
          col.border <- border.allele
        } else if (am) {
          col <- col.filtered
          col.border <- border.filtered
        }
        graphics::polygon(x, y, col = col, border = col.border, lwd = 0.5)
      }
    }
  }
}

#' Plot Sequence Alignments
#'
#' Plot an MSA alignment object.
#'
#' @param alignment MSA alignment object as produced by
#'   \code{\link{align_alleles}}, or character vector of the corresponding
#'   sequences.
#' @param labels custom labels to draw for each entry in \code{alignment}.  By
#'   default it's assumed that \code{align_alleles} was called with
#'   \code{derep=TRUE} and sequences are labeled by number of occurrences.
#' @param include.blanks should blank sequences present in the alignment be
#'   included in the plot?  \code{FALSE} by default.  If TRUE and \code{labels}
#'   is left at the default, the extra axis labels will add up to a full count
#'   of the number of alleles observed.
#' @param ... additional arguments passed to \code{\link[dnaplotr]{plotDNA}}.
#'
#' @return list of the sequence, sequence group, and label character vectors
#'   used in the plot.
#'
#' @seealso \code{\link{align_alleles}}
#'
#' @export
plot_alignment <- function(alignment, labels=NULL, include.blanks=FALSE, ...) {
  # Convert to character and remove blanks if specified
  if (is.character(alignment))
    seqs <- alignment
  else
    seqs <- as.character(alignment)
  if (! include.blanks)
    seqs <- seqs[grep("^-+$", seqs, invert = TRUE)]
  # Create grouping factor using sequence length (just strip out the gap
  # character to get the original length back).  Make sure we're using a
  # consistent order for the sequences, grouping factor, and labels.
  lengths <- nchar(gsub("-", "", seqs))
  ord <- order(lengths)
  seqs <- seqs[ord]
  lengths <- lengths[ord]
  # Add a label for every unique sequence, using the names of the supplied
  # sequences
  if (missing(labels))
    labels <- sapply(strsplit(names(seqs), "_"), "[", 2)
  else
    labels <- labels[ord]
  groups <- paste("  ", lengths, "bp")
  groups <- factor(groups, levels = unique(groups))
  # Make plot
  graphics::par(mar = c(5, 5, 4, 5))
  dnaplotr::plotDNA(seqs,
                    groups = groups,
                    ...)
  # Add faint lines between all sequences
  for (i in seq_along(seqs))
    graphics::abline(h = i + 0.5, col = grDevices::rgb(0.5, 0.5, 0.5, 0.5))
  if (!is.null(labels))
    graphics::axis(4,
         at = 1:length(seqs),
         labels = labels,
         tick = F,
         padj = -2.5,
         cex.axis = 0.6)
  return(list(seqs = seqs, groups = groups, labels = labels))
}

# Distance Matrices -------------------------------------------------------


# A skewed 0 -> 1 scale for color-coding distance tables
make.dist_scale <- function(n) {
  ( (0:n) / n) ** (1 / 3)
}

#' Plot Distance Matrix
#'
#' Plot a heatmap of a distance matrix.
#'
#' @param dist_mat distance matrix as produced by
#'   \code{\link{summarize_dataset}} via \code{\link{make_dist_mat}}.
#' @param num.alleles the maximum number of matching/mis-matching alleles.  Used
#'   to determine color scaling.  Defaults to the highest observed distance in
#'   the matrix.
#' @param dist.display_thresh distance value at or below which distances will be
#'   explicitly drawn on the heatmap.  Above this value only the color-coding
#'   will signify distance.  Use \code{NA} to always show numbers.
#' @param ... additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @seealso \code{\link{make_dist_mat}}
#'
#' @export
plot_dist_mat <- function(dist_mat, num.alleles=max(dist_mat),
                          dist.display_thresh=round(num.alleles * 2 / 3),
                          ...) {
  labels <- matrix(character(length(dist_mat)), nrow = nrow(dist_mat))
  if (is.na(dist.display_thresh))
    dist.display_thresh <- max(dist_mat)
  idx <- dist_mat <= dist.display_thresh
  labels[idx] <- dist_mat[idx]
  diag(labels) <- ""
  dist_scale <- make.dist_scale(num.alleles)
  color <- grDevices::rgb(red = 1, green = dist_scale, blue = dist_scale)

  # Scale font size automatically between min and max values
  fontsize <- min(16, max(4, 17 - 0.11 * max(dim(dist_mat))))

  args <- list(mat = dist_mat,
               color = color,
               display_numbers = labels,
               treeheight_row = 0,
               breaks = 0:num.alleles,
               fontsize = fontsize)
  if (nrow(dist_mat) == ncol(dist_mat)) {
    args <- c(args,
              clustering_distance_rows = stats::as.dist(dist_mat),
              clustering_distance_cols = stats::as.dist(dist_mat))
  } else {
    args <- c(args,
              cluster_cols = FALSE,
              cluster_rows = FALSE)
  }

  do.call(pheatmap::pheatmap, c(args, ...))
}


# Heatmaps ----------------------------------------------------------------

#' Render heatmap of STR attribute across samples and loci
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}} and the name of a column (e.g., Stutter,
#' Homozygous, ProminentSequences), plot a heatmap of the values for that
#' attribute, with sample identifiers on rows and loci on columns.  The
#' attribute will be coerced to numeric.
#'
#' @param results combined results list
#' @param attribute character name of column in results$summary to use.
#' @param label.by vector of column names to use when writing the genotype
#'   summary values on top of the heatmap cells.  Defaults to allele sequence
#'   lengths.
#' @param color vector of colors passed to \code{\link[pheatmap]{pheatmap}}.
#' @param breaks vector of breakpoints passed to
#'   \code{\link[pheatmap]{pheatmap}}, autocalculated by default.
#' @param ... additional arguments to \code{\link[pheatmap]{pheatmap}}.
#'
#' @export
plot_heatmap <- function(results,
                         attribute,
                         label.by = c("Allele1Length", "Allele2Length"),
                         color=c("white", "pink"),
                         breaks=NA,
                         ...) {
  tbl <- summarize_attribute(results$summary, attribute)
  data <- tbl[, - (1:2)] + 0
  tbl.labels <- summarize_genotypes(results$summary, vars = label.by)
  labels <- tbl.labels[, - (1:2)]
  data[is.na(labels)] <- NA
  labels[is.na(labels)] <- ""

  # Handle edge cases where all the values are the same and/or all NA
  if (all(is.na(data)))
    data[, ] <- 0
  if (min(data, na.rm = T) == max(data, na.rm = T))
    breaks <- range(c(0, max(data), 1))

  pheatmap::pheatmap(data,
                     cluster_rows = F,
                     cluster_cols = F,
                     display_numbers = labels,
                     breaks = breaks,
                     color = color,
                     ...)
}

#' Plot heatmap of suspected PCR stutter
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}}, plot a heatmap showing which samples had
#' alleles ignored due to suspected PCR stutter, with sample identifiers on rows
#' and loci on columns.
#'
#' @param results combined results list
#' @param ... additional arguments passed to \code{\link{plot_heatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_heatmap_stutter <- function(results, ...) {
  plot_heatmap(results,
               "Stutter",
               legend = FALSE,
               ...)
}

#' Plot heatmap of suspected PCR artifacts
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}}, plot a heatmap showing which samples had
#' alleles ignored due to a suspected PCR artifact, with sample identifiers on
#' rows and loci on columns.
#'
#' @param results combined results list
#' @param ... additional arguments passed to \code{\link{plot_heatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_heatmap_artifacts <- function(results, ...) {
  plot_heatmap(results,
               "Artifact",
               legend = FALSE,
               ...)
}

#' Plot heatmap of homozygous samples
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}}, plot a heatmap showing which samples appear
#' homozygous, with sample identifiers on rows and loci on columns.
#'
#' @param results combined results list
#' @param ... additional arguments passed to \code{\link{plot_heatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_heatmap_homozygous <- function(results, ...) {
  plot_heatmap(results,
               "Homozygous",
               legend = FALSE,
               ...)
}

#' Plot heatmap of samples with multiple prominent sequences
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}}, plot a heatmap showing samples with more than
#' two prominent sequences in their analysis output, with sample identifiers on
#' rows and loci on columns.
#'
#' @param results combined results list
#' @param ... additional arguments passed to \code{\link{plot_heatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_heatmap_prominent_seqs <- function(results, ...) {
  # Create a color ramp going from white for 0, 1, or 2 prominent seqs, and
  # shades of red for higher numbers.
  color_func <- grDevices::colorRampPalette(c("white", "red"))
  ps <- results$summary[!is.na(results$summary$Allele1Seq), "ProminentSeqs"]
  if (length(ps) == 0) {
    ps <- 0
  }
  # Deep red will only be used if somehow there are a whole lot of extra
  # sequences (say, 8); otherwise it should just go up to a pink color.
  colors <- color_func(max(8, max(ps) + 1))
  # Stay white for 0 - 2
  colors[1:3] <- rep(colors[1], 3)
  # Truncate to actual number of peaks
  colors <- colors[1:max(ps) + 1]
  plot_heatmap(results,
               "ProminentSeqs",
               color = colors,
               breaks = 0:(length(colors)),
               ...)
}

#' Plot heatmap of proportion of allele sequence counts
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}}, plot a heatmap showing the proportion of
#' matching sequences for the identified alleles versus all matching sequences,
#' with sample identifiers on rows and loci on columns.
#'
#' @param results combined results list
#' @param ... additional arguments passed to \code{\link{plot_heatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_heatmap_proportions <- function(results, ...) {
  cts <- results$summary[, c("Allele1Count", "Allele2Count")]
  prop.counted <- rowSums(cts, na.rm = T) / results$summary$CountLocus
  results$summary$ProportionCounted <- prop.counted
  # A color scale going from red at 0 to white at 1, but values skewed toward
  # white.
  color_func <- grDevices::colorRampPalette(c("red", "white"))
  breaks <- seq(0, 1, 0.001) ^ 2.5
  colors <- color_func(length(breaks) - 1)
  plot_heatmap(results,
               "ProportionCounted",
               color = colors,
               breaks = breaks,
               ...)
}


# Counts per Locus Heatmap ------------------------------------------------

#' Plot heatmap of read counts matching each locus primer
#'
#' Given a data frame as produced by \code{\link{tally_cts_per_locus}}, plot a
#' heatmap showing the number of reads matching the forward primer of each locus
#' across samples.  Samples are shown on rows with the reads categorized by
#' locus across columns.
#'
#' @param cts_per_locus data frame as produced by
#'   \code{\link{tally_cts_per_locus}}.
#' @param idx.row Optional vector of sample row indices to use.  (Using this
#'   argument rather than filtering the input allows the same plot scale to be
#'   used across plots.)
#' @param ... additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_cts_per_locus <- function(cts_per_locus, idx.row=NULL, ...) {
  # Switch to log scale
  cts_per_locus[cts_per_locus == 0] <- NA
  cts_per_locus <- log10(cts_per_locus)
  # Break on powers of ten (since we already log10'd above)
  breaks <- 0:ceiling(max(cts_per_locus, na.rm = T))
  color <- viridis::viridis(max(breaks))

  if (! missing(idx.row)) {
    cts_per_locus <- cts_per_locus[idx.row, ]
  }
  pheatmap::pheatmap(cts_per_locus,
                     cluster_rows = F,
                     cluster_cols = F,
                     gaps_col = c(1, 2),
                     color = color,
                     breaks = breaks,
                     legend_breaks = breaks,
                     legend_labels = paste0("10^", breaks),
                     ...)
}
