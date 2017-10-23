# Create summary representations across samples in a dataset.

# between-sample distance matrix
# alignments across identified alleles per locus
# cross-locus contamination reporting
# labels for known alleles? ("180a," "180b," etc.)
# identification of individuals?

# Genotype and Attribute Tables -------------------------------------------

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
  rownames(tbl) <- make.rownames(tbl)
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


# Distances ---------------------------------------------------------------


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
  attr(dist.mat, "Labels") <- rownames(tbl)
  attr(dist.mat, "Size") <- nrow(tbl)
  attr(dist.mat, "Diag") <- attr(dist.mat, "Upper") <- FALSE
  dist.mat <- as.matrix(dist.mat)
  # set the diagonal explicitly, since the distance function may possibly return
  # non-zero values for a distance between an entry and itself.
  diag(dist.mat) <- apply(tbl, 1, function(row) {
    calc.genotype.distance(row[-(1:2)], row[-(1:2)])
  })
  # This all ends up being pretty roundabout.  Should I instead use outer() or
  # something to build the matrix?
  dist.mat
}


# TODO: work with either allele direction
# for test: test symmetry
calc.genotype.distance <- function(g1, g2, na.reject = TRUE) {
  # Valid matches are either at the index in g2 or one less, to match either
  # allele.  So for each position in g2 minus its index that will be -1 or 0.
  #m <- if (na.reject) match(g1, g2, incomparables = NA) else match(g1, g2)
  #sum(! (m - seq_along(g2)) %in% -1:0)
  if (na.reject) {
    g1[is.na(g1)] <- -1
    g2[is.na(g2)] <- -2
  }
  sum(g1 != g2, na.rm = T)
}


# Alignments --------------------------------------------------------------


# create a list of alignments across observed alleles, with one alignment per
# locus and one sequence per allele.
# TODO support dereplication, either here or separately.
align.alleles <- function(dataset.summary) {
  chunks <- split(dataset.summary, dataset.summary$Locus)
  lapply(chunks, function(chunk) {
    alleles <- chunk[, c("Allele1Seq", "Allele2Seq")]
    a1 <- as.character(alleles[,1])
    a2 <- as.character(alleles[,2])
    a2 <- ifelse(is.na(a2), a1, a2)
    names(a1) <- paste(rownames(alleles), 1, sep="_")
    names(a2) <- paste(rownames(alleles), 2, sep="_")
    a <- c(a1, a2)
    a[is.na(a)] <- '-'
    # TODO make this safer.  if msa(...) crashes sink() won't get called.
    sink('/dev/null')
    alignments <- msa::msaClustalW(a,
                                   type="dna",
                                   substitutionMatrix = 'clustalw')
    sink()
    alignments
  })
}
