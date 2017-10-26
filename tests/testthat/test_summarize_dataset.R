prepare_for_summary <- function() {
  data.dir <- tempfile()
  write_seqs(seqs, data.dir)
  dataset <- prepare.dataset(data.dir, '()(\\d+)-([A-Za-z0-9]+).fasta')
  results <- analyze.dataset(dataset, locus_attrs)
  lapply(dataset$Filename, file.remove)
  file.remove(data.dir)
  return(list(dataset=dataset, results=results))
}

results_summary_data <- prepare_for_summary()

# test_summarize.dataset --------------------------------------------------

test_that("summarize.dataset produces additional summaries", {
  # Just a general check of the overall process and binding these results
  # together in a list; more rigorous checks are below since those are public
  # functions too.
  with(results_summary_data, {
    results_mod <- summarize.dataset(results)
    expect_equal(names(results_mod),
                 c(names(results), c("alignments", "dist_mat")))
  })
})

# test_make.dist_mat ------------------------------------------------------

test_that("make.dist_mat produces a valid distance matrix", {
  with(results_summary_data, {
    dist_mat <- make.dist_mat(results$summary)
    # row and column names should both be equal to the sample identifiers
    samps <- levels(dataset$Sample)
    expect_equal(rownames(dist_mat), samps)
    expect_equal(colnames(dist_mat), samps)
    # For our test data samples are far apart from each other, but with no
    # missing entries the diagonal is zero.
    dists <- matrix(7, nrow=length(samps), ncol=length(samps))
    dists[1,3] <- 8
    dists[3,1] <- 8
    diag(dists) <- 0
    expect_true(all(dist_mat == dists))
  })
})

test_that("make.dist_mat_known produces sample-to-individual distance matrix", {
  fail("test not yet implemented")
})

# test_calc.genotype.distance ---------------------------------------------

test_that("calc.genotype.distance scores genotypes correctly", {
  g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
  g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
  d <- calc.genotype.distance(g1, g2)
  expect_equal(d, 4)
})

test_that("calc.genotype.distance handles multiple data types", {
  # As written it should work with character, numeric, whatever-- it just uses
  # the != operator.
  g.num <- c(1,2,3,4)
  g.alpha <- c('a', 'b', 'c', 'd')
  d <- calc.genotype.distance(g.num, g.alpha)
  d.num <- calc.genotype.distance(g.num, g.num)
  d.alpha <- calc.genotype.distance(g.alpha, g.alpha)
  expect_equal(d, 4)
  expect_equal(d.num, 0)
  expect_equal(d.alpha, 0)
})

test_that("calc.genotype.distance is commutative", {
  g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
  g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
  d.fwd <- calc.genotype.distance(g1, g2)
  d.rev <- calc.genotype.distance(g2, g1)
  expect_equal(d.fwd, 4)
  expect_equal(d.rev, 4)
})

test_that("calc.genotype.distance handles NA entries", {
  g1 <- c(5,   5, 200, 204, 37, 39)
  g2 <- c(NA, NA, 204, 208, 37, 37)
  g3 <- c(NA, NA, 200, 204, 37, 39)

  # The two NA entries will by default count towards the distance whether the
  # other genotype has an NA there or not.  If na.reject is FALSE, NAs are
  # treated like any other entry.
  d12 <- calc.genotype.distance(g1, g2)
  d12.na <- calc.genotype.distance(g1, g2, na.reject = FALSE)
  d23 <- calc.genotype.distance(g2, g3)
  d23.na <- calc.genotype.distance(g2, g3, na.reject = FALSE)
  expect_equal(d12,    4)
  expect_equal(d12.na, 2)
  expect_equal(d23,    4)
  expect_equal(d23.na, 2)
})

test_that("calc.genotype.distance invalid genotype lengths", {
  # These are valid but not the same length
  g1 <- c(5,   5, 200, 204)
  g2 <- c(NA, NA, 204, 208, 37, 39)
  # This has an odd length
  g3 <- c(NA, NA, 204, 208, 37)

  expect_warning(calc.genotype.distance(g1, g2))
  expect_warning(calc.genotype.distance(g3, g3))
})


# test_align.alleles ------------------------------------------------------

test_that("align.alleles produces sequence alignments", {
  with(results_summary_data, {
    alignments <- align.alleles(results$summary)
    expect_equal(names(alignments), levels(results$summary$Locus))
    expect_equal(names(as.character(alignments[["A"]])),
                 c("182_1", "194_1", "178_1", "174_2", "162_1"))
    expect_equal(as.character(alignments[["A"]])[[1]],
                 paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGA",
                        "TAGATAGATAGATAGATAGATAGATAGATA",
                        "GATAGATAGATAGATAGATAGATAGATAGA",
                        "TAGATAGATAGATAGATAGATAGATAGATA",
                        "GATAGATAGATAGATAGATAGATAGATAGA",
                        "TAGATAG------------ATAGACACAGT",
                        "TGTGTGAGCCAGTC"))
  })
})

test_that("align.alleles produces per-allele alignments", {
  # By default the allele sequences are dereplicated and then aligned, but this
  # tests the other option.
  with(results_summary_data, {
    alignments <- align.alleles(results$summary, derep = FALSE)
    expect_equal(names(alignments), levels(results$summary$Locus))
    n1 <- paste(rownames(subset(results$summary, Locus == "A")), 1, sep='_')
    n2 <- paste(rownames(subset(results$summary, Locus == "A")), 2, sep='_')
    expect_equal(sort(names(as.character(alignments[["A"]]))),
                 sort(c(n1, n2)))
  })
})
