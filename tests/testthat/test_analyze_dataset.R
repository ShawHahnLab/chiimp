testrds <- function(fname) readRDS(test_path("data", "analyze_dataset", fname))


# Currently if I try to run analyze_dataset with the parallel package the tests
# fail when run under R CMD check but work fine when run interactively.  For now
# I'm just sticking with ncores = 1 here and in test_summarize_dataset.R to
# avoid calling parallel:: functions.
# Possibly relevant:
#  * https://github.com/r-lib/testthat/issues/602
#  * https://github.com/hadley/devtools/issues/1526
#  * https://github.com/r-lib/testthat/issues/86


# test analyze_dataset ----------------------------------------------------


test_that("analyze_dataset produces outputs across samples", {
  dataset <- testrds("dataset.rds")
  seqs <- testrds("seqs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  results_expected <- testrds("results.rds")
  within_tmpdir({
    write_seqs(seqs, "data")
    results <- analyze_dataset(
      dataset, locus_attrs, analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500), ncores = 1)
  })
  expect_equal(results, results_expected)
})

test_that("analyze_dataset names known alleles", {
  # If we gave names for some known allele sequences, do they show up
  # appropriately in the summary table?
  dataset <- testrds("dataset.rds")
  seqs <- testrds("seqs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  known_alleles <- testrds("known_alleles.rds")
  results_expected <- testrds("results_known_alleles.rds")
  within_tmpdir({
    write_seqs(seqs, "data")
    results <- analyze_dataset(
      dataset, locus_attrs, analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500), ncores = 1,
      known_alleles = known_alleles)
  })
  expect_equal(results, results_expected)

  # Check that the resulting allele names match all the expected values
  with(results$summary, {
    expect_equal(Allele1Name, c("different_name_format",
                                "different_name_format",
                                "182-d679e1", "252-27c5bf", "236-321c79",
                                "236-321c79", "280-a", "284-2b3fab",
                                "280-a", "250-5dacee", "266-2aa675",
                                "342-2e88c0"))
    expect_equal(Allele2Name, c("194-fc013a", "178-d84dc0", NA,
                                "216-c0f11a", "240-2a344f", "220-fb9a92",
                                "260-X", "256-c18a06", "276-ea279a",
                                "318-35b7b6", NA, "238-6cc8ff"))
  })

  # Check that the resulting sequence names in the sample data frames match
  # up. We have the one or two called alleles per sample checked here, with
  # one cross-appearance checked below.
  lapply(rownames(results$summary), function(nm) {
    # First called allele for these cases should always be the first seq in
    # each table.
    expect_equal(results$summary[nm, "Allele1Name"],
                 results$samples[[nm]]$SeqName[1])
    # Second called allele, if present, will be below the first somewhere.
    # Remaining seqs will be unnamed.
    if (! results$summary[nm, "Homozygous"]) {
      idx <- match(results$summary[nm, "Allele2Seq"],
                   results$samples[[nm]]$Seq)
      expect_equal(results$summary[nm, "Allele2Name"],
                   results$samples[[nm]]$SeqName[idx])
    }
  })

  # One particular case: 3-B showed a stutter-rejected sequence that's the
  # called allele for another sample.
  expect_equal(results$samples[["3-B"]]$SeqName[2], "220-fb9a92")
})

test_that("analyze_dataset handles missing loci", {
  # If there are locus names in dataset$Locus that are not present in the
  # rownames of locus_attrs, it should throw an error.
  dataset <- testrds("dataset.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  # the names are case-sensitive!
  dataset$Locus[dataset$Locus == "A"] <- "a"
  dataset$Locus[dataset$Locus == "B"] <- "b"
  expect_error(
    analyze_dataset(
      dataset, locus_attrs, analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500), ncores = 1),
    "ERROR: Locus names in dataset not in attributes table: a, b")
})

test_that("analyze_dataset warns of empty input files", {
  # If we have no reads at all right from the start, we should warn the user.
  dataset <- testrds("dataset.rds")
  seqs <- testrds("seqs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  within_tmpdir({
    write_seqs(seqs, "data")
    # empty out one file
    fps <- list.files("data", full.names = TRUE)
    unlink(fps[1])
    touch(fps[1])
    msg <- capture.output({
      results <- analyze_dataset(
        dataset, locus_attrs, analysis_opts = list(fraction.min = 0.05),
        summary_opts = list(counts.min = 500), ncores = 1)
    }, type = "message")
    msg_exp <- "WARNING: Zero reads for 1 of 12 data files"
    expect_true(length(grep(msg_exp, msg)) == 1)
  })
})
