# Create summary representations across samples in a dataset.

# Full Dataset Summary ----------------------------------------------------

#' Add further summaries to analyzed dataset
#'
#' Take a results list as produced by \code{\link{analyze_dataset}} and add
#' additional entries for inter-sample and inter-locus analyses.
#'
#' @details
#' Additional entries in the returned list:
#'   * alignments: inter-allele alignments for each locus, from
#'     \code{\link{align_alleles}}.
#'   * dist_mat: inter-sample distance matrix, from \code{\link{make_dist_mat}}.
#'   * dist_mat_known: if genotypes.known is given, this distance matrix of
#'     sample-to-individual values will be present, from
#'     \code{\link{make_dist_mat_known}}.
#' @md
#'
#' @param results list containing summary data frame and sample-specific data
#'   frames as produced by \code{\link{analyze_dataset}}.
#' @param genotypes.known optional data frame of known genotypes that should be
#'   compared to the observed genotypes in the results, as loaded by
#'   \code{\link{load_genotypes}}.  If provided \code{dist_mat_known} will be
#'   present in the output.
#'
#' @return expanded list with additional summaries.
#'
#' @export
summarize_dataset <- function(results, genotypes.known=NULL) {
  results$alignments <- align_alleles(results$summary)
  results$dist_mat <- make_dist_mat(results$summary)
  if (!missing(genotypes.known) & !is.null(genotypes.known)) {
    results$dist_mat_known <- make_dist_mat_known(results$summary,
                                                  genotypes.known)
    results$genotypes.known <- genotypes.known
  }
  return(results)
}

# Distances ---------------------------------------------------------------

# TODO this is all super roundabout.  rewrite to use the outer() approach for
# both and just have one version of the bulk of the method.

#' Make distance matrix for set of samples
#'
#' Compare the genotype of each sample to each other sample in a given results
#' summary data frame, and create a between-sample distance matrix.  This will
#' be a symmetric matrix, unlike what \code{\link{make_dist_mat_known}}
#' produces.
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param dist.func function to calculate inter-sample distances.  Should take
#'   two vectors, one for each genotype, with two values per locus corresponding
#'   to the two alleles.
#'
#' @return matrix of cross-sample distance values
#'
#' @export
make_dist_mat <- function(results_summary,
                          dist.func=calc_genotype_distance) {
  tbl <- summarize_genotypes(results_summary)
  # The entire vector covers all combinations of rows in tbl, filling in a
  # triangle of what would be a distance matrix.
  distances <- utils::combn(nrow(tbl), 2,
        function(nr) {
          dist.func(tbl[nr[1], - (1:2)],
                    tbl[nr[2], - (1:2)])
        })
  # Trick from SO to funnel it through a dist object to get it into matrix form.
  # https://stackoverflow.com/a/5598824/6073858
  dist.mat <- distances
  class(dist.mat) <- "dist"
  attr(dist.mat, "Labels") <- rownames(tbl)
  attr(dist.mat, "Size") <- nrow(tbl)
  attr(dist.mat, "Diag") <- attr(dist.mat, "Upper") <- FALSE
  dist.mat <- as.matrix(dist.mat)
  # set the diagonal explicitly, since the distance function may possibly return
  # non-zero values for a distance between an entry and itself.
  diag(dist.mat) <- apply(tbl, 1, function(row) {
    dist.func(row[- (1:2)], row[- (1:2)])
  })
  # This all ends up being pretty roundabout.  Should I instead use outer() or
  # something to build the matrix?
  dist.mat
}

#' Make distance matrix for samples to known genotypes
#'
#' Compare the genotype of each sample in a given results summary to all known
#' genotypes in the supplied table, and create a sample-to-known-genotype
#' distance matrix.
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param genotypes.known data frame of known genotypes, with one row per
#'   individual per locus, and columns "Name", "Locus", "Allele1Seq",
#'   "Allele2Seq".  See \code{\link{load_genotypes}}.
#' @param dist.func function to calculate inter-sample distances.  Should take
#'   two vectors, one for each genotype, with two values per locus corresponding
#'   to the two alleles.
#'
#' @return matrix of sample-to-individual distance values, with individuals from
#'   genotypes.known on rows and samples on columns.
#'
#' @export
make_dist_mat_known <- function(results_summary,
                                genotypes.known,
                                dist.func=calc_genotype_distance) {
  tbl <- summarize_genotypes(results_summary)
  tbl.known <- summarize_genotypes_known(genotypes.known, tbl)
  distances <- outer(rownames(tbl),
                     rownames(tbl.known),
                     function(nr, nr.known) {
                       apply(cbind(nr, nr.known), 1,
                             function(row) {
                               dist.func(tbl[row[1], - (1:2)],
                                         tbl.known[row[2], - (1:2)])
                             })
  })
  rownames(distances) <- rownames(tbl)
  colnames(distances) <- rownames(tbl.known)
  distances
}

#' Calculate distance between two genotypes
#'
#' Given two vectors representing genotypes, produce a scalar value for the
#' distance between them.  The vectors should be sorted the same such that the
#' each two values of each vector represent two alleles of the same locus.  In
#' this implementation the distance is simply the number of mismatches between
#' the alleles for each locus, accounting for either order of alleles.
#'
#' @param g1 vector for genotype.
#' @param g2 vector for genotype.
#' @param na.reject logical; should NA entries be treated as mismatches?  TRUE
#'   by default.
#'
#' @return numeric distance score for the pair of input genotypes.
#'
#' @export
calc_genotype_distance <- function(g1, g2, na.reject = TRUE) {
  g1 <- unlist(g1)
  g2 <- unlist(g2)
  if (length(g1) %% 2 != 0 || length(g2) %% 2 != 0) {
    warning("Odd length for input genotype; truncating.")
    g1 <- g1[1:(length(g1) - length(g1) %% 2)]
    g2 <- g2[1:(length(g2) - length(g2) %% 2)]
  }
  if (length(g1) != length(g2)) {
    warning("Input genotype length mismatch; truncating.")
    g1 <- g1[1:min(length(g1), length(g2))]
    g2 <- g2[1:min(length(g1), length(g2))]
  }

  if (na.reject) {
    g1[is.na(g1)] <- -1
    g2[is.na(g2)] <- -2
  }
  alleles1 <- matrix(g1, ncol = 2, byrow = TRUE)
  alleles2 <- matrix(g2, ncol = 2, byrow = TRUE)
  alleles <- cbind(alleles1, alleles2)
  alleles
  sum(apply(alleles, 1,
            function(row) min(sum(row[1:2] != row[3:4], na.rm = T),
                              sum(row[2:1] != row[3:4]), na.rm = T)))
}

#' Find closest matches in distance matrix
#'
#' Given a distance matrix with samples on rows and names on columns, return a
#' list of vectors with the closest-matching names for each sample.
#'
#' @param dist_mat matrix of distance values, such as produced by
#'   \code{\link{make_dist_mat}} and \code{\link{make_dist_mat_known}}.
#' @param range optional numeric for distances to each set of nearby names,
#'   relative to the closest match.
#' @param maximum optional numeric maximum value for any distance.
#'
#' @return list of named vectors containing distances for each sample.
#'
#' @export
find_closest_matches <- function(dist_mat, range=2, maximum=8) {
  entries <- lapply(1:nrow(dist_mat), function(nr) {
    m <- min(dist_mat[nr, ])
    nearby <- dist_mat[nr, dist_mat[nr, ] < m + range &
                         dist_mat[nr, ] < maximum, drop = F]
    nearby <- nearby[1, order(nearby), drop = F]
    nm <- colnames(nearby)
    nearby <- nearby[1, ]
    names(nearby) <- nm
    nearby
  })
  names(entries) <- rownames(dist_mat)
  entries
}


# Alignments --------------------------------------------------------------

#' Align allele sequences across loci
#'
#' Create a list of alignments across observed alleles, with one alignment per
#' locus present in the given summary, and by default one entry for each unique
#' sequence for each locus.
#'
#' @param results_summary cross-sample summary data frame.
#' @param derep logical; should allele sequences be dereplicated prior to
#'   alignment?  TRUE by default.
#' @param ... additional arguments passed to \code{\link[msa]{msa}}.
#'
#' @return list of MSA alignment objects, one per locus.
#'
#' @export
align_alleles <- function(results_summary, derep=TRUE, ...) {
  chunks <- split(results_summary, droplevels(results_summary$Locus))
  lapply(chunks, function(chunk) {
    alleles <- chunk[, c("Allele1Seq", "Allele2Seq")]
    a1 <- as.character(alleles[, 1])
    a2 <- as.character(alleles[, 2])
    a2 <- ifelse(is.na(a2), a1, a2)
    names(a1) <- paste(rownames(alleles), 1, sep = "_")
    names(a2) <- paste(rownames(alleles), 2, sep = "_")
    a <- c(a1, a2)
    # If there are no sequences, skip the alignment and record NA for this
    # locus.  (The msa function doesn't do this sort of checking itself,
    # apparently.
    if (all(is.na(a)))
      return(NULL)
    # Turn NAs into empty strings
    a[is.na(a)] <- ""
    # Dereplicate identical sequences, if specified
    if (derep) {
      tbl <- table(a)
      n <- unname(tbl)
      a <- names(tbl)
      names(a) <- paste(nchar(a), n, sep = "_")
    }
    # If there are any other missing sequences, put a stub in place so msa still
    # runs without complaints.
    a[a == ""] <- "-"
    # msa() generates a bunch of text on standard output and I can't see any
    # options to turn that off.  Using a workaround here.
    tryCatch({
      sink(fp_devnull)
      alignments <- msa::msaClustalW(a,
                                     type = "dna",
                                     substitutionMatrix = "clustalw",
                                     ...)
    },
    finally = {
      sink()
    })
    alignments
  })
}

# Genotype and Attribute Tables -------------------------------------------

#' Create Summary Genotype Table
#'
#' Tabulate genotypes across an analyzed dataset, with samples on rows and loci
#' on columns.  Which attributes will actually be used for the table is
#' configurable.
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param vars vector of column names in \code{results_summary} to use for the
#'   summary.  Defaults to the sequence content for the two alleles.
#'
#' @return data frame of genotypes across samples and loci.
#'
#' @export
summarize_genotypes <- function(results_summary,
                                vars=c("Allele1Seq", "Allele2Seq")) {
  combo <- results_summary[, c("Sample", "Replicate", "Locus", vars)]
  combo[, vars[1]] <- as.character(combo[, vars[1]])
  combo[, vars[2]] <- as.character(combo[, vars[2]])
  # Repeat alleles for homozygous cases
  combo[, vars[2]] <- ifelse(results_summary$Homozygous,
                             combo[, vars[1]],
                             combo[, vars[2]])
  # Re-order to put shorter alleles first, if we've been called with vars for
  # integers (otherwise we're just sorting text)
  if (is.character(combo[[vars[1]]])) {
    idx <- nchar(combo[[vars[1]]]) > nchar(combo[[vars[2]]]) |
      ( nchar(combo[[vars[1]]]) == nchar(combo[[vars[2]]]) &
          combo[[vars[1]]] > combo[[vars[2]]] )
  } else {
    idx <- combo[[vars[1]]] > combo[[vars[2]]]
  }
  idx[is.na(idx)] <- F
  combo[idx, vars] <- combo[idx, rev(vars)]
  # Mark those cases that are not homozygous but just happen to have the same
  # short name (i.e. allele sequence length)
  combo[, vars[2]] <- ifelse( !(results_summary$Homozygous) &
                                combo[, vars[1]] == combo[, vars[2]],
                              paste0(combo[, vars[2]], "*"),
                              combo[, vars[2]])
  # Reshape into wide table
  combo$ID <- paste(combo$Sample, combo$Replicate, sep = "-")
  tbl <- stats::reshape(combo, v.names = vars, idvar = "ID",
                        timevar = "Locus", direction = "wide")
  tbl <- tbl[, -3]
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each = 2),
                       c(1, 2),
                       sep = "_")
  colnames(tbl) <- c("Sample", "Replicate", allele_cols)
  tbl <- tbl[order_entries(tbl), ]
  rownames(tbl) <- make_rownames(tbl)
  tbl
}

#' Create Summary Genotype Table for Known Genotypes
#'
#' Tabulate genotypes across an analyzed dataset, with samples on rows and loci
#' on columns.  Which attributes will actually be used for the table is
#' configurable.
#'
#' @param genotypes_known table of known individuals' genotypes as loaded via
#'   \code{\link{load_genotypes}}.
#' @param tbl_genotypes existing genotype summary table to use for locus
#'   selection and column ordering.
#'
#' @return data frame of genotypes across individuals and loci.
#'
#' @export
summarize_genotypes_known <- function(genotypes_known, tbl_genotypes=NULL) {
  # Kludgy workaround to make summarize_genotypes handle a different sort of
  # data frame
  genotypes_known$Replicate <- NA
  genotypes_known$Sample <- genotypes_known$Name
  genotypes_known <- genotypes_known[, -match("Name",
                                              colnames(genotypes_known))]
  genotypes_known$Homozygous <- is.na(as.character(genotypes_known$Allele2Seq))
  tbl_known <- summarize_genotypes(genotypes_known)
  if (!missing(tbl_genotypes))
    tbl_known <- tbl_known[, colnames(tbl_genotypes)]
  tbl_known
}

#' Create Summary Attribute Table
#'
#' Tabulate a single arbitrary attribute across loci, assuming repeats by two
#' for the alleles.  This is used for color-coding summary heatmaps (see
#' \code{\link{plot_heatmap}}) on top of the attribute values, like Homozgyous
#' or ProminentSeqs.
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param attrib name of column from \code{results_summary} to tabulate.
#' @param repeats number of times to repeat each column, for matching with two
#'   alleles per entry.
#'
#' @return data frame of attribute across samples and loci.
#'
#' @export
summarize_attribute <- function(results_summary, attrib, repeats = 2) {
  attrib_rep <- rep(attrib, repeats)
  combo <- results_summary[, c("Sample", "Replicate", "Locus", attrib_rep)]
  combo$ID <- paste(combo$Sample, combo$Replicate, sep = "-")
  # reshape into wide table
  tbl <- stats::reshape(combo,
                        v.names = make.names(attrib_rep, unique = T),
                        idvar = "ID",
                        timevar = "Locus", direction = "wide")
  tbl <- tbl[, -3]
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each = repeats),
                       1:repeats,
                       sep = "_")
  colnames(tbl) <- c("Sample", "Replicate", allele_cols)
  tbl <- tbl[order_entries(tbl), ]
  rownames(tbl) <- make_rownames(tbl)
  tbl
}

# Counts per Locus --------------------------------------------------------

#' Create Counts-per-Locus Table
#'
#' Produce a table of read counts matching (by primer) each locus per sample.
#'
#' @param results list containing summary data frame and sample-specific data
#'   frames as produced by \code{\link{analyze_dataset}}.
#'
#' @return Data frame of read counts with samples on rows and loci on columns.
#'   Addition columns at the left specify the number of reads matching any known
#'   locus (\code{Total}) and the number of reads matching the expected locus
#'   (\code{Matching}).
#'
#' @export
tally_cts_per_locus <- function(results) {
  # Create table of counts of sequences that match each possible locus across
  # samples.  Only include loci we expect from the metadata, rather than any
  # known in locus_attrs.
  .levels <- match(rownames(results$locus_attrs), results$summary$Locus)
  .levels <- results$summary$Locus[.levels[!is.na(.levels)]]
  tbl <- do.call(rbind, lapply(results$data, function(d) {
    d$MatchingLocus <- factor(d$MatchingLocus,
                              levels = .levels)
    sapply(split(d$Count, d$MatchingLocus), sum)
  }))
  # Make some extra columns for total sequences and sequences matching the
  # expected locus.  Bind these to the original data to force the heatmap to use
  # a uniform scale.
  cols.match <- results$summary[rownames(tbl), "Locus"]
  tbl.anno <- data.frame(Total = rowSums(tbl),
                         Matching = sapply(seq_along(cols.match),
                            function(i) tbl[i, as.character(cols.match[i])]))
  tbl <- cbind(tbl.anno, tbl)
  tbl
}

# Locus Performance -------------------------------------------------------

# TODO: functions to evaluate the level of genotyping success and allele
# diversity across loci.
