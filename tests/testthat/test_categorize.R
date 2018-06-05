context("Test categorization of genotyping results for known samples")

with(test_data, {


  # test match_known_genotypes --------------------------------------------


  test_that("match_known_genotypes arranges genotypes to match samples", {
    # A basic test to make sure we get the right data frame back and it contains
    # the expected values for a simple case.
    with(results_summary_data, {
      # Leaving the third sample unidentified.
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      output <- match_known_genotypes(results$summary, genotypes_known)
      expect_equal(colnames(output),
                            c("CorrectAllele1Seq", "CorrectAllele2Seq"))
      # This is a simplistic test; every sample's genotype will match correct
      # genotype, and with the same ordering of the alleles.  We'll just make
      # sure that there are at least some non-NA entries and that for those the
      # alleles match.
      expect_true(any(! is.na(output$CorrectAllele1Seq)))
      expect_true(all(is.na(output$CorrectAllele1Seq) |
            results$summary$Allele1Seq == output$CorrectAllele1Seq))
      expect_true(any(! is.na(output$CorrectAllele2Seq)))
      expect_true(all(is.na(output$CorrectAllele2Seq) |
                        results$summary$Allele2Seq == output$CorrectAllele2Seq))
      # The way we set up the identities, the first sample should be the second
      # known individual (skip four loci, go to 5th entry).  For example
      expect_equal(unname(unlist(
        output[1, c("CorrectAllele1Seq", "CorrectAllele2Seq")])),
        unname(unlist(
          genotypes_known[5, c("Allele1Seq", "Allele2Seq")])))
      # What about that third sample, that didn't match any known individual?
      # There should be NA for all of those entries.
      expect_true(all(is.na(output[results$summary$Sample == "3", ])))
    })
  })

  test_that("match_known_genotypes orders alleles to match samples", {
    # match_known_genotypes is supposed to flip the order of the alleles given
    # in the known genotypes data frame if that would match up to the order
    # given in the analysis results.  By default we've matched them up in the
    # correct order already here.
    with(results_summary_data, {
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      # A simple flip for the first case.
      results$summary[1, c("Allele1Seq", "Allele2Seq")] <-
        results$summary[1, c("Allele2Seq", "Allele1Seq")]
      # For the second, screw up the allele sequences.
      results$summary[2, "Allele1Seq"] <- "ACGT"
      results$summary[2, "Allele2Seq"] <- "GTCA"
      # For the third, we'll go to another ID002 sample and drop an allele.  But
      # let's drop the *first* as ordered in the known genotypes table.
      results$summary[4, "Allele1Seq"] <- results$summary[4, "Allele2Seq"]
      results$summary[4, "Allele2Seq"] <- NA
      output <- match_known_genotypes(results$summary, genotypes_known)
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
  })


  # test categorize_genotype_results --------------------------------------


  test_that("categorize_genotype_results categorizes results", {
    # A basic test to make sure we get the right factor back with the expected
    # values for a simple case.
    with(results_summary_data, {
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      output <- match_known_genotypes(results$summary, genotypes_known)
      results$summary <- cbind(results$summary, output)
      categories <- categorize_genotype_results(results$summary)
      # The output will use four defined categories, in this order.
      expect_equal(levels(categories),
                   c("Correct", "Dropped Allele", "Blank", "Incorrect"))
      # In this simple case we have some Correct, and some with no comment since
      # we don't have an individual in our set with that identifier.
      expect_equal(as.character(categories),
                   rep(c("Correct", "Correct", NA), 4))
    })
  })

  test_that("categorize_genotype_results identifies dropped alleles", {
    with(results_summary_data, {
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      # We'll drop one allele each from the first two entries, trying both the
      # first and second alleles.
      results$summary[1, "Allele2Seq"] <- NA
      results$summary[2, "Allele1Seq"] <- NA
      output <- match_known_genotypes(results$summary, genotypes_known)
      results$summary <- cbind(results$summary, output)
      categories <- categorize_genotype_results(results$summary)
      expect_equal(as.character(categories[1:2]),
                   rep("Dropped Allele", 2))
    })
  })

  test_that("categorize_genotype_results identifies blanks", {
    with(results_summary_data, {
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      # We'll drop one full allele pair from the first entry, and also the
      # third.  The third has no known individual to match to, so that result
      # should still stay NA, though.
      results$summary[1, c("Allele1Seq", "Allele2Seq")] <- NA
      results$summary[3, c("Allele1Seq", "Allele2Seq")] <- NA
      output <- match_known_genotypes(results$summary, genotypes_known)
      results$summary <- cbind(results$summary, output)
      categories <- categorize_genotype_results(results$summary)
      expect_equal(as.character(categories[c(1, 3)]),
                   c("Blank", NA))
    })
  })

  test_that("categorize_genotype_results identifies incorrect alleles", {
    with(results_summary_data, {
      results$summary$Name <- c("ID002", "ID001")[
        as.integer(results$summary$Sample)]
      # We'll change the allele sequence for one allele of the first entry, one
      # of the second, and both of the fourth.  All of these should then be
      # Incorrect.
      results$summary[1, "Allele1Seq"] <- "ACTG"
      results$summary[2, "Allele2Seq"] <- "ACTG"
      results$summary[4, "Allele2Seq"] <- "ACTG"
      output <- match_known_genotypes(results$summary, genotypes_known)
      results$summary <- cbind(results$summary, output)
      categories <- categorize_genotype_results(results$summary)
      expect_equal(as.character(categories[1:4]),
                   c("Incorrect", "Incorrect", NA, "Incorrect"))
    })
  })

})
