# Create summary representations across samples in a dataset.

# between-sample distance matrix
# alignments across identified alleles per locus
# cross-locus contamination reporting
# labels for known alleles? ("180a," "180b," etc.)
# identification of individuals?

# Tabulate genotypes across loci.  Which attributes will actually be used for
# the table is configurable.
summarize.genotypes <- function(dataset.summary,
                                vars=c('Allele1Length', 'Allele2Length')) {
  combo <- dataset.summary[, c('Sample', 'Replicate', 'Locus', vars)]
  # Repeat alleles for homozygous cases
  combo[, vars[2]] <- ifelse(dataset.summary$Homozygous,
                             combo[, vars[1]],
                             combo[, vars[2]])
  # Re-order to put shorter alleles first (TODO generalize this beyond integer
  # comparison)
  idx <- combo[[vars[1]]] > combo[[vars[2]]]
  idx[is.na(idx)] <- F
  combo[idx, vars] <- combo[idx, rev(vars)]
  # Mark those cases that are not homozygous but just happen to have the same
  # short name (i.e. allele sequence length)
  combo[, vars[2]] <- ifelse( !(dataset.summary$Homozygous) &
                                dataset.summary[, vars[1]] == dataset.summary[, vars[2]],
                              paste0(combo[, vars[2]], '*'),
                              combo[, vars[2]])
  # Reshape into wide table
  tbl <- reshape(combo, v.names = vars, idvar = 'Sample',
                 timevar = 'Locus', direction = 'wide')
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each=2),
                       c(1,2),
                       sep = '_')
  colnames(tbl) <- c('Sample', 'Replicate', allele_cols)
  tbl
}

# tabulate an arbitrary attribute across loci, assuming repeats by two for the
# alleles.  (use this for color-coding summary heatmaps on top of the attribute
# values, like Homozgyous or ProminentSeqs.)
summarize.attribute <- function(dataset.summary, attrib, repeats = 2) {
  attrib <- rep(attrib, repeats)
  combo <- dataset.summary[, c('Sample', 'Replicate', 'Locus', attrib)]
  # reshape into wide table
  tbl <- reshape(combo, v.names = colnames(attrib_data), idvar = 'Sample',
                 timevar = 'Locus', direction = 'wide')
  allele_cols <- paste(rep(as.character(unique(combo$Locus)), each = repeats),
                       1:repeats,
                       sep = '_')
  colnames(tbl) <- c('Sample', 'Replicate', allele_cols)
  tbl
}

# Create a distance matrix for every combination of individuals in the
# supplied dataset results summary.
make.dist_mat <- function(dataset.summary,
                          dist.func=calc.genotype.distance) {
  tbl <- summarize.genotypes(dataset.summary)
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
  attr(dist.mat, "Labels") <- tbl$Sample
  attr(dist.mat, "Size") <- nrow(tbl)
  attr(dist.mat, "Diag") <- attr(dist.mat, "Upper") <- FALSE
  dist.mat <- as.matrix(dist.mat)
  # set the diagonal explicitly, since the distance function may possibly return
  # non-zero values for a distance between an entry and itself.
  diag(dist.mat) <- apply(tbl, 1, function(row) {
    calc.genotype.distance(row[-(1:2)], row[-(1:2)])
  })
  dist.mat
}

# TODO update this to consider re-arranged alleles.
calc.genotype.distance <- function(g1, g2, na.reject = TRUE) {
  if (na.reject) {
    g1[is.na(g1)] <- -1
    g2[is.na(g2)] <- -2
  }
  sum(g1 != g2, na.rm = T)
}
