# Create summary representations across samples in a dataset.

# between-sample distance matrix
# alignments across identified alleles per locus
# cross-locus contamination reporting
# labels for known alleles? ("180a," "180b," etc.)
# identification of individuals?


# Full Dataset Summary ----------------------------------------------------

#' Add further summaries to analyzed dataset
#'
#' Take a results list as produced by \code{analyze.dataset} and add additional
#' entries for inter-sample and inter-locus analyses:
#'
#' @details
#' Additional entries in the returned list:
#'   * alignments: inter-allele alignments for each locus
#'   * dist_mat: inter-sample distance matrix
#'   * dist_mat_known: if genotypes.known is given, this distance matrix of
#'     sample-to-individual values will be present.
#' @md
#'
#' @param results list containing summary data frame and sample-specific data
#'   frames as produced by \code{analyze.dataset}.
#' @param genotypes.known optional data frame of known genotypes that should be
#'   compared to the observed genotypes in the results.  If provided
#'   dist_mat_known will be present in the output.
#'
#' @return expanded list with additional summaries.
#'
#' @export
summarize.dataset <- function(results, genotypes.known=NULL) {
  results$alignments <- align.alleles(results$summary)
  results$dist_mat <- make.dist_mat(results$summary)
  if (!missing(genotypes.known)) {
    results$dist_mat_known <- make.dist_mat_known(results$summary,
                                                  genotypes.known)
  }
  return(results)
}

# Genotype and Attribute Tables -------------------------------------------

# Tabulate genotypes across loci.  Which attributes will actually be used for
# the table is configurable.
summarize.genotypes <- function(results_summary,
                                vars=c('Allele1Seq', 'Allele2Seq')) {
  combo <- results_summary[, c('Sample', 'Replicate', 'Locus', vars)]
  # Repeat alleles for homozygous cases
  if (class(combo[, vars[1]]) == 'factor')
    combo[, vars[1]] <- as.character(combo[, vars[1]])
  if (class(combo[, vars[2]]) == 'factor')
    combo[, vars[2]] <- as.character(combo[, vars[2]])
  combo[, vars[2]] <- ifelse(results_summary$Homozygous,
                             combo[, vars[1]],
                             combo[, vars[2]])
  # Re-order to put shorter alleles first, if we've been called with vars for
  # integers (otherwise we're just sorting text)
  if (is.character(combo[[vars[1]]])) {
    idx <- nchar(combo[[vars[1]]]) > nchar(combo[[vars[2]]])
  } else {
    idx <- combo[[vars[1]]] > combo[[vars[2]]]
  }
  idx[is.na(idx)] <- F
  combo[idx, vars] <- combo[idx, rev(vars)]
  # Mark those cases that are not homozygous but just happen to have the same
  # short name (i.e. allele sequence length)
  combo[, vars[2]] <- ifelse( !(results_summary$Homozygous) &
                                combo[, vars[1]] == combo[, vars[2]],
                              paste0(combo[, vars[2]], '*'),
                              combo[, vars[2]])
  # Reshape into wide table
  tbl <- reshape(combo, v.names = vars, idvar = 'Sample',
                 timevar = 'Locus', direction = 'wide')
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each=2),
                       c(1,2),
                       sep = '_')
  colnames(tbl) <- c('Sample', 'Replicate', allele_cols)
  rownames(tbl) <- make.rownames(tbl)
  tbl
}

# tabulate an arbitrary attribute across loci, assuming repeats by two for the
# alleles.  (use this for color-coding summary heatmaps on top of the attribute
# values, like Homozgyous or ProminentSeqs.)
summarize.attribute <- function(results_summary, attrib, repeats = 2) {
  attrib_rep <- rep(attrib, repeats)
  combo <- results_summary[, c('Sample', 'Replicate', 'Locus', attrib_rep)]
  combo$ID <- paste(combo$Sample, combo$Replicate, sep='-')
  # reshape into wide table
  tbl <- reshape(combo,
                 v.names = make.names(attrib_rep, unique = T),
                 idvar = 'ID',
                 timevar = 'Locus', direction = 'wide')
  tbl <- tbl[, -3]
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each = repeats),
                       1:repeats,
                       sep = '_')
  colnames(tbl) <- c('Sample', 'Replicate', allele_cols)
  rownames(tbl) <- make.rownames(tbl)
  tbl
}


# Distances ---------------------------------------------------------------

# TODO this is all super roundabout.  rewrite to use the outer() approach for
# both and just have one version of the bulk of the method.

#' Make distance matrix for set of samples
#'
#' Compare the genotype of each sample to each other sample in a given results
#' summary data frame, and create a between-sample distance matrix.  This will
#' be a symmetric matrix, unlike what \code{make.dist_mat_known} produces.
#'
#' @param results_summary cross-sample summary data frame
#' @param dist_func function to calculate inter-sample distances.  Should take
#'   two vectors, one for each genotype, with two values per locus corresponding
#'   to the two alleles.
#'
#' @return matrix of cross-sample distance values
#'
#' @export
make.dist_mat <- function(results_summary,
                          dist.func=calc.genotype.distance) {
  tbl <- summarize.genotypes(results_summary)
  # The entire vector covers all combinations of rows in tbl, filling in a
  # triangle of what would be a distance matrix.
  distances <- combn(nrow(tbl), 2,
        function(nr) {
          dist.func(tbl[nr[1], -(1:2)],
                    tbl[nr[2], -(1:2)])
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
    dist.func(row[-(1:2)], row[-(1:2)])
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
#' @param results_summary cross-sample summary data frame
#' @param genotypes.known data frame of known genotypes, with one row per
#'   individual per locus, and columns "Name", "Locus", "Allele1Seq",
#'   "Allele2Seq".
#' @param dist_func function to calculate inter-sample distances.  Should take
#'   two vectors, one for each genotype, with two values per locus corresponding
#'   to the two alleles.
#'
#' @return matrix of sample-to-individual distance values, with individuals from
#'   genotypes.known on rows and samples on columns.
#'
#' @export
make.dist_mat_known <- function(results_summary,
                                genotypes.known,
                                dist.func=calc.genotype.distance) {
  tbl <- summarize.genotypes(results_summary)
  genotypes.known$Replicate <- NA
  genotypes.known$Sample <- genotypes.known$Name
  genotypes.known$Homozygous <- as.character(genotypes.known$Allele1Seq) ==
    as.character(genotypes.known$Allele2Seq)
  tbl.known <- summarize.genotypes(genotypes.known)
  distances <- outer(rownames(tbl),
                     rownames(tbl.known),
                     function(nr, nr.known) {
                       apply(cbind(nr, nr.known), 1,
                             function(row) {
                               dist.func(tbl[row[1], -(1:2)],
                                         tbl.known[row[2], -(1:2)])
                             })
  })
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
calc.genotype.distance <- function(g1, g2, na.reject = TRUE) {
  g1 <- unlist(g1)
  g2 <- unlist(g2)
  if (length(g1)%%2 != 0 || length(g2)%%2 != 0) {
    warning("Odd length for input genotype; truncating.")
    g1 <- g1[1:(length(g1)-length(g1)%%2)]
    g2 <- g2[1:(length(g2)-length(g2)%%2)]
  }
  if (length(g1) != length(g2)) {
    warning("Input genotype length mismatch; truncating.")
    g1 <- g1[1:min(length(g1), length(g2))]
    g2 <- g2[1:min(length(g1), length(g2))]
  }

  # Valid matches are either at the index in g2 or one less, to match either
  # allele.  So for each position in g2 minus its index that will be -1 or 0.
  #m <- if (na.reject) match(g1, g2, incomparables = NA) else match(g1, g2)
  #sum(! (m - seq_along(g2)) %in% -1:0)
  if (na.reject) {
    g1[is.na(g1)] <- -1
    g2[is.na(g2)] <- -2
  }
  alleles1 <- matrix(g1, ncol = 2, byrow = TRUE)
  alleles2 <- matrix(g2, ncol = 2, byrow = TRUE)
  alleles <- cbind(alleles1, alleles2)
  alleles
  # TODO this mis-counts things like c(1,1) versus c(1,2).  Maybe just go back
  # to the forward/reverse check from before.
  sum(apply(alleles, 1, function(row) min(sum(row[1:2] != row[3:4], na.rm=T),
                                          sum(row[2:1] != row[3:4]), na.rm=T) ))
  #sum(g1 != g2, na.rm = T)
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
#' @param ... additional arguments passed to \code{msa}.
#'
#' @return list of MSA alignment objects, one per locus.
#'
#' @export
align.alleles <- function(results_summary, derep=TRUE, ...) {
  chunks <- split(results_summary, results_summary$Locus)
  lapply(chunks, function(chunk) {
    alleles <- chunk[, c("Allele1Seq", "Allele2Seq")]
    a1 <- as.character(alleles[,1])
    a2 <- as.character(alleles[,2])
    a2 <- ifelse(is.na(a2), a1, a2)
    names(a1) <- paste(rownames(alleles), 1, sep="_")
    names(a2) <- paste(rownames(alleles), 2, sep="_")
    a <- c(a1, a2)
    a[is.na(a)] <- '-'
    if (derep) {
      tbl <- table(a)
      n <- unname(tbl)
      a <- names(tbl)
      names(a) <- paste(nchar(a), n,sep='_')
    }

    # TODO make this safer.  if msa(...) crashes sink() won't get called.  Use
    # tryCatch?
    sink('/dev/null')
    alignments <- msa::msaClustalW(a,
                                   type="dna",
                                   substitutionMatrix = 'clustalw',
                                   ...)
    sink()
    alignments
  })
}


# Locus Performance -------------------------------------------------------

# TODO: functions to evaluate the level of genotyping success and allele
# diversity across loci.
