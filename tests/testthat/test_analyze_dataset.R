# test analyze_dataset ----------------------------------------------------

# Cases to cover: fastq, fasta, fasta.gz, fastq.gz, multiplex, ...?
# At this point we should have all the basics covered, and can cover the further
# analysis in summarize_dataset.

# Currently if I try to run analyze_dataset with the parallel package the tests
# fail when run under R CMD check but work fine when run interactively.  For now
# I'm just sticking with ncores = 1 here and in test_summarize_dataset.R to
# avoid calling parallel:: functions.
# Possibly relevant:
#  * https://github.com/r-lib/testthat/issues/602
#  * https://github.com/hadley/devtools/issues/1526
#  * https://github.com/r-lib/testthat/issues/86

test_that("analyze_dataset processes samples correctly", {
  # The general case for analyze_dataset.
  data.dir <- tempfile()
  write_seqs(seqs, data.dir)
  # prepare_dataset tested separately in test_io.R
  dataset <- prepare_dataset(data.dir, '()(\\d+)-([A-Za-z0-9]+).fasta')
  results <- analyze_dataset(dataset, locus_attrs,
                             summary_args = list(
                               fraction.min = 0.05,
                               counts.min = 500), nrepeats = 3, ncores = 1)
  with(results, {
    # These should just match what we fed in via dataset.
    expect_equal(summary$Filename,  dataset$Filename)
    expect_equal(summary$Replicate, dataset$Replicate)
    expect_equal(summary$Sample,    dataset$Sample)
    expect_equal(summary$Locus,     dataset$Locus)
    # Skipping allele seqs here, but we should already be checking more than
    # enough and have that checked in test_summarize_sample.
    expect_equal(summary$Allele1Count, c(2794, 2041, 2909, 3035, 2902, 2288,
                                         2803, 2394, 3971, 2489, 2656, 4082))
    expect_equal(summary$Allele1Length, c(272, 276, 276, 278, 346, 318,
                                          162, 182, 174, 244, 212, 224))
    expect_equal(summary$Allele2Count, c(1055, 1402, 1004, 1142, 1088, 1270,
                                         1300, 2003,   NA, 1304, 1059,   NA))
    expect_equal(summary$Allele2Length, c(260, 288, 252, 342, 314, 350,
                                          194, 178,  NA, 220, 220,  NA))
    h <- logical(12)
    h[c(9, 12)] <- TRUE
    expect_equal(summary$Homozygous, h)
    s <- logical(12)
    s[c(9, 12)] <- TRUE
    expect_equal(summary$Stutter, s)
    expect_equal(summary$CountTotal, integer(12)+5000)
    expect_equal(summary$CountLocus, integer(12)+4500)
    expect_equal(summary$ProminentSeqs,  c(4, 3, 4, 2, 3, 3, 3, 2, 1, 3, 4, 1))
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)
  })

})

test_that("analyze_dataset runs single-threaded", {
  # analyze_dataset should run properly with ncores=1 (just slower).
  skip("test not yet implemented")
})

test_that("analyze_dataset warns of missing loci", {
  # If there are locus names in dataset$Locus that are not present in the
  # rownames of locus_attrs, it should throw a warning.
  skip("test not yet implemented")
})
