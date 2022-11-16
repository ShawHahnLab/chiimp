testrds <- function(fname) readRDS(test_path("data", "categorize", fname))

# test match_known_genotypes ----------------------------------------------


test_that("match_known_genotypes arranges genotypes to match samples", {
  # A basic test to make sure we get the right data frame back and it contains
  # the expected values for a simple case.
  # Here the results summary table has a Name column with two samples with
  # names assigned and the third without.  We should get CorrectAllele1Seq and
  # CorrectAllele2Seq entries for samples 1 and 2 but not 3.
  results_summary <- testrds("results_summary.rds")
  genotypes_known <- testrds("genotypes_known.rds")
  output_expected <- testrds("match_known_genotypes.rds")
  output <- match_known_genotypes(results_summary, genotypes_known)
  expect_equal(output, output_expected)
})

test_that("match_known_genotypes orders alleles to match samples", {
  # match_known_genotypes is supposed to flip the order of the alleles given
  # in the known genotypes data frame if that would match up to the order
  # given in the analysis results.  By default we've matched them up in the
  # correct order already here.
  results_summary <- testrds("results_summary.rds")
  genotypes_known <- testrds("genotypes_known.rds")
  results_summary$Name <- c("ID002", "ID001")[
    as.integer(results_summary$Sample)]
  # A simple flip for the first case.
  results_summary[1, c("Allele1Seq", "Allele2Seq")] <-
    results_summary[1, c("Allele2Seq", "Allele1Seq")]
  # For the second, screw up the allele sequences.
  results_summary[2, "Allele1Seq"] <- "ACGT"
  results_summary[2, "Allele2Seq"] <- "GTCA"
  # For the third, we'll go to another ID002 sample and drop an allele.  But
  # let's drop the *first* as ordered in the known genotypes table.
  results_summary[4, "Allele1Seq"] <- results_summary[4, "Allele2Seq"]
  results_summary[4, "Allele2Seq"] <- NA
  output <- match_known_genotypes(results_summary, genotypes_known)
  # Now, the first one should match in *reverse* order to the expected
  # individual.
  expect_equal(unname(unlist(
                output[1, c("CorrectAllele1Seq", "CorrectAllele2Seq")])),
               unname(unlist(
                 genotypes_known[5, c("Allele2Seq", "Allele1Seq")])))
  # Second case, the ordering should be left as-is since neither matched.
  # This is the first individual in the table, this time.
  expect_equal(unname(unlist(
    output[2, c("CorrectAllele1Seq", "CorrectAllele2Seq")])),
    unname(unlist(
      genotypes_known[1, c("Allele1Seq", "Allele2Seq")])))
  # Third case, dropped allele; the order should be flipped again.  back to
  # IND002, second locus.
  expect_equal(unname(unlist(
    output[4, c("CorrectAllele1Seq", "CorrectAllele2Seq")])),
    unname(unlist(
      genotypes_known[6, c("Allele2Seq", "Allele1Seq")])))
})


# test categorize_genotype_results ----------------------------------------


test_that("categorize_genotype_results categorizes results", {
  # A basic test to make sure we get the right factor back with the expected
  # values for a simple case.
  results_summary <- testrds("results_summary_matched.rds")
  categories_expected <- testrds("categories.rds")
  categories <- categorize_genotype_results(results_summary)
  # The output should use four defined categories, ordered as "Correct",
  # "Dropped Allele", "Blank", and "Incorrect".  In this case some should be
  # labeled Correct and some NA (since we have a sample without a matched
  # individual)
  expect_equal(categories, categories_expected)
})

test_that("categorize_genotype_results identifies dropped alleles", {
  # This example is missing one allele each for the first two entries, the first
  # allele for the first and the second allele for the second.
  results_summary <- testrds("results_summary_matched_drop.rds")
  categories_expected <- testrds("categories_drop.rds")
  categories <- categorize_genotype_results(results_summary)
  # This time we should have two cases of "Dropped Allele" in the vector.
  expect_equal(categories, categories_expected)
})

test_that("categorize_genotype_results identifies blanks", {
  # This example is missing a whole allele pair for the first and third entries.
  # The third has no know individual to match to, so that result should still
  # stay NA, though.
  results_summary <- testrds("results_summary_matched_drop.rds")
  categories_expected <- testrds("categories_drop.rds")
  categories <- categorize_genotype_results(results_summary)
  # Here we should see one Blank in the vector
  expect_equal(categories, categories_expected)
})

test_that("categorize_genotype_results identifies incorrect alleles", {
  results_summary <- testrds("results_summary.rds")
  genotypes_known <- testrds("genotypes_known.rds")
  results_summary$Name <- c("ID002", "ID001")[
    as.integer(results_summary$Sample)]
  # We'll change the allele sequence for one allele of the first entry, one
  # of the second, and both of the fourth.  All of these should then be
  # Incorrect.
  results_summary[1, "Allele1Seq"] <- "ACTG"
  results_summary[2, "Allele2Seq"] <- "ACTG"
  results_summary[4, "Allele2Seq"] <- "ACTG"
  output <- match_known_genotypes(results_summary, genotypes_known)
  results_summary <- cbind(results_summary, output)
  categories <- categorize_genotype_results(results_summary)
  expect_equal(as.character(categories[1:4]),
               c("Incorrect", "Incorrect", NA, "Incorrect"))
})
