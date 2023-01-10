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
    unname(pair)
  }
  )), stringsAsFactors = FALSE)
}

#' Wide table of allele names vs loci
#'
#' Allele pairs are shown in a standardized order with homozygous entries shown
#' twice.  NA allele names in the input are converted to empty strings while NA
#' is given for missing sample/locus combinations.
#'
#' @param data data frame containing Allele1Name and Allele2Name columns such as
#'   the first list item produced by \code{\link{analyze_dataset}}.  If allele
#'   names are not yet present call \code{\link{name_alleles_in_table}}.
#' @param extra_cols names or index values of additional columns from input data
#'   frame to be kept in output data frame.  These should be consistent across
#'   loci for a given entry.
#'
#' @return wide format data frame with sample entries on rows and loci on
#'   columns.  An ID column will label sample entries by whichever columns were
#'   provided in the input (see \code{\link{make_entry_id}}).  Empty strings are
#'   given for NA allele names, while NA is given for any implicitly missing
#'   locus/sample combinations.
#'
#' @export
tabulate_allele_names <- function(data, extra_cols = NULL) {
  if (length(extra_cols)) {
    badcols <- extra_cols[! extra_cols %in% colnames(data)]
    if (length(badcols)) {
      stop(
        paste("undefined extra columns in tabulate_allele_names: ",
              badcols, collapse = " "))
    }
  }
  # Order and replicate (for homozygous) the allele names
  nms <- normalize_alleles(data[, c("Allele1Name", "Allele2Name")])
  # Create unique (aside from Locus) identifiers for each entry
  id <- make_entry_id(data[, -match("Locus", colnames(data))])
  # Our normalized and ordered long-format data frame to be reshaped.
  long <- data.frame(
    ID = id,
    Locus = data$Locus,
    data[, extra_cols, drop = FALSE],
    nms,
    stringsAsFactors = FALSE)
  long <- long[order_entries(long), ]
  # Switch to wide format, putting the allele names per locus across columns
  # (along with ID and whatever extra_cols were given).
  tbl <- stats::reshape(
    long, v.names = c("V1", "V2"), idvar = "ID", timevar = "Locus",
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
  # If extra columns were given, re-order output rows using those
  if (length(extra_cols)) {
    tbl <- tbl[order_entries(tbl[, extra_cols, drop = FALSE]), ]
  }
  tbl
}

#' Create genotype summary table
#'
#' Report the genotypes present in a processed dataset in a concise data frame.
#' This will arrange the allele names into a wide-format table with unique
#' samples on rows and loci on columns, do some automatic cleanup on the
#' columns, and show closest-matching individuals per entry, if given.  All NA
#' entries are replaced with blank strings or optionally (for NA Replicates or
#' untested sample/locus combinations) other custom placeholder text.
#'
#' @param results list of results data as produced by \code{analyze_dataset}.
#' @param na.replicates text to replace NA entries with for the Replicates
#'     column.
#' @param na.alleles text to replace NA entries with for the allele names
#' @param closest list of closest matches as produced by
#'   \code{\link{find_closest_matches}}.
#'
#' @return data frame showing summary of genotypes.
#' @export
report_genotypes <- function(
    results, na.replicates = "", na.alleles = "", closest = NULL) {
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

  # If we have no replicates drop that column.  Otherwise put placeholder text
  # for any NA replicate entries.
  if (all(is.na(tbl$Replicate)))
    tbl <- tbl[, -2]
  else
    tbl$Replicate[is.na(tbl$Replicate)] <- na.replicates

  # Put placeholder text for any untested sample/locus combinations
  # (This is a clumsy way of handling different columns differently, and is
  # probably a hint that more logic handled in the long-format data frames would
  # be better, but this can be a stopgap before some reorganization at some
  # point.)
  locus_cols <- do.call(
    paste0,
    expand.grid(unique(results$summary$Locus), c("_1", "_2")))
  for (colnm in colnames(tbl)) {
      if (colnm %in% locus_cols) {
        tbl[[colnm]][is.na(tbl[[colnm]])] <- na.alleles
      }
  }

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
#' @return data frame showing summary of sample genotypes with interleaved
#'   genotypes for similar known individuals.
#'
#' @export
report_idents <- function(results, closest, na.replicates = "") {
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
plot_alignment <- function(
    alignment, labels = NULL, include.blanks = FALSE, ...) {
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
         at = seq_len(length(seqs)),
         labels = labels,
         tick = FALSE,
         padj = -2.5,
         cex.axis = 0.6)
  return(list(seqs = seqs, groups = groups, labels = labels))
}

# Distance Matrices -------------------------------------------------------


# A skewed 0 -> 1 scale for color-coding distance tables
make.dist_scale <- function(n) {
  ((0:n) / n) ** (1 / 3)
}

#' Plot Distance Matrix
#'
#' Plot a heatmap of a distance matrix.
#'
#' @param dist_mat distance matrix as produced by
#'   \code{\link{summarize_dataset}} via \code{\link{make_dist_mat}}.
#' @param num.alleles the maximum number of matching/mismatching alleles.  Used
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
plot_dist_mat <- function(
    dist_mat, num.alleles = max(dist_mat),
    dist.display_thresh = round(num.alleles * 2 / 3), ...) {
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

  args <- list(
    mat = dist_mat,
    color = color,
    display_numbers = labels,
    treeheight_row = 0,
    breaks = 0:num.alleles,
    fontsize = fontsize)
  if (nrow(dist_mat) == ncol(dist_mat)) {
    args <- c(args,
              list(clustering_distance_rows = stats::as.dist(dist_mat)),
              list(clustering_distance_cols = stats::as.dist(dist_mat)))
  } else {
    args <- c(args,
              cluster_cols = FALSE,
              cluster_rows = FALSE)
  }

  do.call(pheatmap::pheatmap, c(args, list(...)))
}


# Heatmaps ----------------------------------------------------------------

#' Render heatmap of STR attribute across samples and loci
#'
#' Given a cross-sample summary data frame as produced by
#' \code{\link{analyze_dataset}} and the name of a column (e.g., \code{Stutter},
#' \code{Homozygous}, \code{ProminentSequences}), plot a heatmap of the values
#' for that attribute, with sample identifiers on rows and loci on columns.  The
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
plot_heatmap <- function(
    results, attribute, label.by = c("Allele1Length", "Allele2Length"),
    color = c("white", "pink"), breaks = NA, ...) {
  tbl <- summarize_attribute(results$summary, attribute)
  data <- tbl[, - (1:2)] + 0
  tbl.labels <- summarize_genotypes(results$summary, vars = label.by)
  labels <- tbl.labels[, - (1:2)]
  data[is.na(labels)] <- NA
  labels[is.na(labels)] <- ""

  # Handle edge cases where all the values are the same and/or all NA
  if (all(is.na(data)))
    data[, ] <- 0
  if (min(data, na.rm = TRUE) == max(data, na.rm = TRUE))
    breaks <- range(c(0, max(data, na.rm = TRUE), 1))

  pheatmap::pheatmap(
    data, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = labels,
    breaks = breaks, color = color, ...)
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
  plot_heatmap(results, "Stutter", legend = FALSE, ...)
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
  plot_heatmap(results, "Artifact", legend = FALSE, ...)
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
  plot_heatmap(results, "Homozygous", legend = FALSE, ...)
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
  prop.counted <- rowSums(cts, na.rm = TRUE) / results$summary$CountLocus
  results$summary$ProportionCounted <- prop.counted
  # A color scale going from red at 0 to white at 1, but values skewed toward
  # white.
  color_func <- grDevices::colorRampPalette(c("red", "white"))
  breaks <- seq(0, 1, 0.001) ^ 2.5
  colors <- color_func(length(breaks) - 1)
  plot_heatmap(
    results, "ProportionCounted", color = colors, breaks = breaks, ...)
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
#' @param render Should the plot be drawn to the display device?
#' @param ... additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
plot_cts_per_locus <- function(
    cts_per_locus, idx.row = NULL, render = TRUE, ...) {
  # Switch to log scale
  cts_per_locus[cts_per_locus == 0] <- NA
  cts_per_locus <- log10(cts_per_locus)
  # Break on powers of ten (since we already log10'd above)
  if (all(is.na(cts_per_locus))) {
    breaks <- 0:1 # handle the all-zero case
  } else {
    breaks <- 0:ceiling(max(cts_per_locus, na.rm = TRUE))
  }
  color <- viridis::viridis(max(breaks))

  if (! missing(idx.row)) {
    cts_per_locus <- cts_per_locus[idx.row, ]
  }
  pheatmap::pheatmap(
    cts_per_locus,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col = c(1, 2),
    color = color,
    breaks = breaks,
    legend_breaks = breaks,
    legend_labels = paste0("10^", breaks),
    silent = ! render,
    ...)
  invisible()
}
