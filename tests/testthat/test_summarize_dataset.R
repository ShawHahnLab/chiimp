testrds <- function(fn) readRDS(test_path("data", "summarize_dataset", fn))


# test_summarize_dataset --------------------------------------------------


test_that("summarize_dataset produces additional summaries", {
  results <- testrds("results.rds")
  results_expected <- testrds("results_mod.rds")
  results_mod <- summarize_dataset(results)
  expect_equal(results_mod, results_expected)
})

test_that("summarize_dataset works with known genotypes", {
  # Several additional things should happen if genotypes_known is given to
  # summarize_dataset.  The extra data frame should be appended to the results
  # list, and a new distance matrix should be present comparing samples to
  # known genotypes.
  # For cases where a Name column is present in the dataset data frame at the
  # start (and so is now present in the results summary data frame) the
  # summary data frame will be expanded with three additional columns.
  # Two (CorrectAllele1Seq and CorrectAllele2Seq) will track the correct
  # sequences, ordered to match the called alleles where possible.
  # GenotypeResult will provide a category for each case (row): Correct,
  # Incorrect, Dropped Allele, and Blank.
  results <- testrds("results_known.rds")
  results_expected <- testrds("results_known_mod.rds")
  genotypes_known <- testrds("genotypes_known.rds")
  results_mod <- summarize_dataset(results, genotypes_known = genotypes_known)
  expect_equal(results_mod, results_expected)
})


# test_make_dist_mat ------------------------------------------------------


test_that("make_dist_mat produces a valid distance matrix", {
  results <- testrds("results.rds")
  dist_mat_expected <- testrds("dist_mat.rds")
  dist_mat <- make_dist_mat(results$summary)
  expect_equal(dist_mat, dist_mat_expected)
  # row and column names should both be equal to the sample identifiers
  samps <- unique(results$summary$Sample)
  expect_equal(rownames(dist_mat), samps)
  expect_equal(colnames(dist_mat), samps)
  # For our test data samples are far apart from each other, but with no
  # missing entries the diagonal is zero.
  dists <- matrix(7, nrow = length(samps), ncol = length(samps))
  diag(dists) <- 0
  expect_true(all(dist_mat == dists))
})

test_that(
  "make_dist_mat_known produces sample-to-individual distance matrix", {
  skip("test not yet implemented")
})


# test_calc_genotype_distance ---------------------------------------------


test_that("calc_genotype_distance scores genotypes correctly", {
  g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
  g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
  d <- calc_genotype_distance(g1, g2)
  expect_equal(d, 4)
})

test_that("calc_genotype_distance handles multiple data types", {
  # As written it should work with character, numeric, whatever-- it just uses
  # the not-equal operator.
  g.num <- c(1, 2, 3, 4)
  g.alpha <- c("a", "b", "c", "d")
  d <- calc_genotype_distance(g.num, g.alpha)
  d.num <- calc_genotype_distance(g.num, g.num)
  d.alpha <- calc_genotype_distance(g.alpha, g.alpha)
  expect_equal(d, 4)
  expect_equal(d.num, 0)
  expect_equal(d.alpha, 0)
})

test_that("calc_genotype_distance is commutative", {
  g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
  g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
  d.fwd <- calc_genotype_distance(g1, g2)
  d.rev <- calc_genotype_distance(g2, g1)
  expect_equal(d.fwd, 4)
  expect_equal(d.rev, 4)
})

test_that("calc_genotype_distance handles NA entries", {
  g1 <- c(5,   5, 200, 204, 37, 39)
  g2 <- c(NA, NA, 204, 208, 37, 37)
  g3 <- c(NA, NA, 200, 204, 37, 39)

  # The two NA entries will by default count towards the distance whether the
  # other genotype has an NA there or not.  If na_reject is FALSE, NAs are
  # treated like any other entry.
  d12 <- calc_genotype_distance(g1, g2)
  d12.na <- calc_genotype_distance(g1, g2, na_reject = FALSE)
  d23 <- calc_genotype_distance(g2, g3)
  d23.na <- calc_genotype_distance(g2, g3, na_reject = FALSE)
  expect_equal(d12,    4)
  expect_equal(d12.na, 2)
  expect_equal(d23,    4)
  expect_equal(d23.na, 2)
})

test_that("calc_genotype_distance invalid genotype lengths", {
  # These are valid but not the same length
  g1 <- c(5,   5, 200, 204)
  g2 <- c(NA, NA, 204, 208, 37, 39)
  # This has an odd length
  g3 <- c(NA, NA, 204, 208, 37)

  expect_warning(calc_genotype_distance(g1, g2))
  expect_warning(calc_genotype_distance(g3, g3))
})


# test_align_alleles ------------------------------------------------------


test_that("align_alleles produces sequence alignments", {
  results <- testrds("results.rds")
  alignments_expected <- testrds("alignments.rds")
  alignments <- align_alleles(results$summary)
  expect_equal(alignments, alignments_expected)
  expect_equal(names(alignments), levels(results$summary$Locus))
  expect_equal(names(as.character(alignments[["A"]])),
               c("182_2", "194_1", "162_2", "178_1"))
  expect_equal(as.character(alignments[["A"]])[[1]],
               paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGA",
                      "TAGATAGATAGATAGATAGATAGATAGATA",
                      "GATAGATAGATAGATAGATAGATAGATAGA",
                      "TAGATAGATAGATAGATAGATAGATAGATA",
                      "GATAGATAGATAGATAGATAGATAGATAGA",
                      "TAGATAG------------ATAGACACAGT",
                      "TGTGTGAGCCAGTC"))
})

test_that("align_alleles produces per-allele alignments", {
  # By default the allele sequences are dereplicated and then aligned, but
  # this tests the other option.
  results <- testrds("results.rds")
  alignments_expected <- testrds("alignments_no_derep.rds")
  alignments <- align_alleles(results$summary, derep = FALSE)
  expect_equal(alignments, alignments_expected)
  expect_equal(names(alignments), levels(results$summary$Locus))
  n1 <- paste(rownames(subset(results$summary, Locus == "A")), 1, sep = "_")
  n2 <- paste(rownames(subset(results$summary, Locus == "A")), 2, sep = "_")
  expect_equal(sort(names(as.character(alignments[["A"]]))),
               sort(c(n1, n2)))
})

test_that("align_alleles works for empty sequences", {
  # Empty sequences should be handled automatically before the msa function is
  # called.
  results <- testrds("results.rds")
  # Empty out one entry for locus A
  idx <- which(results$summary$Locus == "A")
  results$summary[idx[1], "Allele1Seq"] <- NA
  results$summary[idx[1], "Allele2Seq"] <- NA
  # Locus A's alignment should be OK, just with one stub entry
  alignments <- align_alleles(results$summary)
  expect_true(any(grepl("^-+$", as.character(alignments$A))))
})

test_that("align_alleles works for empty sequence sets", {
  # Completely empty sequence lists for a given locus are a special case; for
  # those we should just get "NA" for the whole entry.
  results <- testrds("results.rds")
  # Empty out all sequences for locus A
  idx <- results$summary$Locus == "A"
  results$summary[idx, "Allele1Seq"] <- NA
  results$summary[idx, "Allele2Seq"] <- NA
  # This should still work, just with NA for A's alignment
  alignments <- align_alleles(results$summary)
  expect_equal(alignments$A, NULL)
})

test_that("align_alleles works for identical sequence sets", {
  # Apparently msa() fails when given a single sequence (derep=TRUE for
  # align_alleles), throwing a "There is an invalid aln file!" error.
  # I'll do its validating for it and account for this inside align_alleles.
  results <- testrds("results.rds")
  # Overwrite Locus A seqs with same content
  idx <- which(results$summary$Locus == "A")
  results$summary[idx, "Allele1Seq"] <- results$summary[idx[1],
                                                        "Allele1Seq"]
  results$summary[idx, "Allele2Seq"] <- NA
  alignments <- align_alleles(results$summary)
  # This should still work, just with the same sequence back for A's
  # alignment
  expect_equal(unname(alignments$A),
               as.character(results$summary[idx[1], "Allele1Seq"]))
})


# test_tally_cts_per_locus ------------------------------------------------


test_that("tally_cts_per_locus counts matches per locus per sample", {
  results <- testrds("results.rds")
  cts_expected <- testrds("cts_per_locus.rds")
  cts_observed <- tally_cts_per_locus(results)
  # We should see 5000 reads total, and then 17 split away for each of the
  # four other loci.  Rownames should match sample row IDs.
  expect_equal(cts_observed, cts_expected)
})

test_that("tally_cts_per_locus uses only loci present in samples", {
  # What if we drop Locus 2 from the results, as though we'd only analyzed
  # three loci even though we have four defined in the attributes table?
  # Locus 2 should not show up in the cts_per_locus data frame.
  # NOTE: Currently the Total column shows the total across the loci shown,
  # not the grand total of reads in the raw file.
  results <- testrds("results.rds")
  cts_expected <- testrds("cts_per_locus.rds")
  # We should see 5000 reads total, and then 17 split away for each of the
  # four other loci.  Rownames should match sample row IDs.
  ids_rm <- c("1-2", "2-2", "3-2")
  cts_expected <- subset(
    cts_expected, ! rownames(cts_expected) %in% ids_rm, select = -`2`)
  cts_expected$Total <- rowSums(cts_expected[, 3:5])
  results$summary <- results$summary[-match(ids_rm,
                                            rownames(results$summary)), ]
  results$samples <- results$samples[-match(ids_rm,
                                            names(results$samples))]
  results$summary$Locus <- droplevels(results$summary$Locus)
  # (This is a slightly goofy way of testing this outcome since each data
  # frame in results$files still will have "2" in the MatchingLocus column,
  # but it should work.)
  cts_observed <- tally_cts_per_locus(results)
  expect_equal(cts_observed, cts_expected)
})
