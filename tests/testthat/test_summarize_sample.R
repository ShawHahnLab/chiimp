testrds <- function(fname) readRDS(test_path("data", "summarize_sample", fname))


# test summarize_sample ---------------------------------------------------

test_that("summarize_sample summarizes sample attributes", {
  sample_data <- testrds("sample_data.rds")
  sample_summary_expected <- testrds("sample_summary.rds")
  sample_summary <- summarize_sample(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample handles completely empty sample data", {
  sample_data <- testrds("sample_data_empty.rds")
  sample_summary_expected <- testrds("sample_summary_empty.rds")
  sample_summary <- summarize_sample(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  # A zero-row data frame in should work here, just producing a list with a
  # bunch of 0/NA/FALSE entries.
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample handles empty sequences in input sample data", {
  sample_data <- testrds("sample_data_stubs.rds")
  sample_summary_expected <- testrds("sample_summary_stubs.rds")
  sample_summary <- summarize_sample(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  # Nothing should change in the output, except that we zeroed out 90 reads
  # that would otherwise get counted (the rest were already set to off-target
  # junk during setup).
  expect_equal(sample_summary, sample_summary_expected)
})


test_that("summarize_sample marks stutter removal", {
  # (different sample here than above)
  sample_data <- testrds("sample_data_stutter.rds")
  sample_summary_expected <- testrds("sample_summary_stutter.rds")
  sample_summary <- summarize_sample(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  # Here an entry that could have been called as an allele was rejected as
  # presumed stutter, so the Stutter value in the summary list should now be
  # TRUE.
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample marks removal of multiple artifact types", {
  # There are separate columns in the output for ambiguous/stutter/artifact
  # sequences and each one can independently be TRUE if a sequence was skipped
  # due to that particular aspect.  (Previously only the first would be
  # marked, and it was incorrectly assumed that a secondary sequence was
  # always the one that may have been filtered.)
  # In this version half of the "correct" reads are split off with an
  # ambiguous sequence.  This also leaves some existing stutter reads in
  # place.
  sample_data <- testrds("sample_data_multi_artifact.rds")
  sample_summary_expected <- testrds("sample_summary_multi_artifact.rds")
  sample_summary <- summarize_sample(sample_data, min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample handles multiple stutter sequences", {
  # If multiple candidate allele sequences are marked as potential stutter,
  # they should all be skipped, not just the first.
  # In this version there are two high-abundance stutter sequences, and both
  # should be recognized as stutter.
  sample_data <- testrds("sample_data_multi_stutter.rds")
  sample_summary_expected <- testrds("sample_summary_multi_stutter.rds")
  sample_summary <- summarize_sample(sample_data, min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})


test_that("summarize_sample handles ambiguous sequences", {
  # If sequences contain unexpected characters these should be interpreted as
  # ambiguous and filtered.
  # In this version the second-highest-count sequence is modified to include an
  # ambiguous base near the end.
  sample_data <- testrds("sample_data_ambig.rds")
  sample_summary_expected <- testrds("sample_summary_ambig.rds")
  sample_summary <- summarize_sample(sample_data, min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample rejects low-count samples", {
  # Here we check that the filtered-counts-thresholding is applied, by checking
  # a case where the counts to a low number.  This should still report some
  # stats but should leave out the allele1/allele2 information
  sample_data <- testrds("sample_data_low.rds")
  sample_summary_expected <- testrds("sample_summary_low.rds")
  sample_summary <- summarize_sample(sample_data, min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})


test_that("summarize_sample counts prominent sequences", {
  sample_datas <- testrds("sample_datas_b.rds")
  sample_summaries_expected <- testrds("sample_summaries_b.rds")
  sample_summaries <- lapply(
    sample_datas, summarize_sample, list(Locus = "B"), 500)
  # Should have 4, 2, and 3 prominent peaks, respectively
  expect_equal(sample_summaries[[1]], sample_summaries_expected[[1]])
  expect_equal(sample_summaries[[2]], sample_summaries_expected[[2]])
  expect_equal(sample_summaries[[3]], sample_summaries_expected[[3]])
})

test_that("summarize_sample works with vector for sample attrs", {
  sample_data <- testrds("sample_data.rds")
  sample_summary_expected <- testrds("sample_summary.rds")
  # sample_attrs argument is a vector this time but the outcome should be
  # exactly the same
  sample_summary <- summarize_sample(
    sample_data, c(Locus = "A"), min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})


# test summarize_sample_guided --------------------------------------------


sample.summary.cols <- c("Allele1Seq", "Allele1Count",
                         "Allele1Length", "Allele2Seq",
                         "Allele2Count", "Allele2Length",
                         "Homozygous", "Ambiguous", "Stutter", "Artifact",
                         "CountTotal", "CountLocus", "ProminentSeqs")

test_that("summarize_sample warns of missing locus name", {
  # summarize_sample() should be able to tell if an invalid (as per the
  # earlier processing) locus name is given, because it won't be in the levels
  # of the MatchingLocus factor of the data frame from analyze_seqs().  I
  # think this would only come about when locus names given by
  # prepare.dataset() don't match what's in locus_attrs.
  skip("test not yet implemented")
})

check_seqs1a_summary <- function(data,
                                 count_locus = 4466,
                                 allele1_count = 2783,
                                 allele2_count = 1290,
                                 ord = 1:2) {
  expect_equal(names(data), sample.summary.cols)
  with(data, {
    alleles <- c(Allele1Seq, Allele2Seq)[ord]
    counts  <- c(Allele1Count, Allele2Count)[ord]
    lengths <- c(Allele1Length, Allele2Length)[ord]
    expect_equal(alleles[1],
                 gsub("[\n ]*", "",
                      "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                      ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                      TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACACAG
                      TTGTGTGAGCCAGTC"))
    expect_equal(alleles[2],
                 gsub("[\n ]*", "",
                      "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                      ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                      TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT
                      AGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC"))
    expect_equal(counts[1],     allele1_count)
    expect_equal(lengths[1],    162)
    expect_equal(counts[2],     allele2_count)
    expect_equal(lengths[2],    194)
    expect_equal(Homozygous,    FALSE)
    expect_equal(Ambiguous,     FALSE)
    expect_equal(Stutter,       FALSE)
    expect_equal(Artifact,      FALSE)
    expect_equal(CountTotal,    5000)
    expect_equal(CountLocus,    count_locus)
    expect_equal(ProminentSeqs, 2)
  })
}

test_that("summarize_sample_guided summarizes sample attributes", {
  # This should have the same behavior as summarize_sample above.
  sample_data <- testrds("sample_data_guided.rds")
  sample_summary_expected <- testrds("sample_summary_guided.rds")
  sample_summary <- summarize_sample_guided(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample_guided handles completely empty sample data", {
  # This should have the same behavior as summarize_sample above.
  sample_data <- testrds("sample_data_empty.rds")
  sample_summary_empty_expected <- testrds(
    "sample_summary_empty.rds")
  sample_summary <- summarize_sample_guided(
    sample_data, list(Locus = "A"), min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_empty_expected)
})

test_that("summarize_sample_guided uses expected_lengths", {
  # This should give allele sequences matching the given expected_lengths,
  # including order.
  # Flipping the order of alleles here to check that aspect
  sample_data <- testrds("sample_data_guided_flip.rds")
  sample_summary_expected <- testrds("sample_summary_guided_flip.rds")
  sample_attrs <- list(
    Locus = "A", ExpectedLength1 = 194, ExpectedLength2 = 162)
  sample_summary <- summarize_sample_guided(sample_data, sample_attrs,
                                            min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample_guided ignores min_locus_reads for expected lengths", {
  # This should give allele sequences matching the given expected_lengths,
  # including order, despite the min_locus_reads value.
  # Flip the order of alleles here to check that aspect
  sample_data <- testrds("sample_data_guided_flip.rds")
  sample_summary_expected <- testrds("sample_summary_guided_flip.rds")
  sample_attrs <- list(
    Locus = "A", ExpectedLength1 = 194, ExpectedLength2 = 162)
  sample_summary <- summarize_sample_guided(sample_data, sample_attrs,
                                            min_locus_reads = 5000)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample_guided ignores min_locus_reads for one exp. length", {
  # This should give allele sequences matching the given expected_lengths,
  # including order, despite the min_locus_reads value.
  sample_data <- testrds("sample_data_guided_one_len.rds")
  sample_summary_expected <- testrds("sample_summary_guided_one_len.rds")
  sample_attrs <- list(
    Locus = "A", ExpectedLength1 = 194, ExpectedLength2 = 194)
  sample_summary <- summarize_sample_guided(
    sample_data, sample_attrs, min_locus_reads = 5000)
  # There should be one and only one called allele.
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample_guided uses min_locus_reads if no expected lengths", {
  # This should not report alleles if total filtered read count is below
  # min_locus_reads threshold.
  sample_attrs <- list(Locus = "A", ExpectedLength1 = NA, ExpectedLength2 = NA)
  sample_data <- testrds("sample_data_guided_no_lens.rds")
  sample_summary_expected <- testrds("sample_summary_guided_no_lens.rds")
  sample_summary <- summarize_sample_guided(
    sample_data, sample_attrs, min_locus_reads = 5000)
  expect_equal(sample_summary, sample_summary_expected)
})

test_that("summarize_sample_guided works with vector for sample attrs", {
  # It's "supposed" to be a list, but when used with functions like apply with
  # a data frame it gets munged into a vector.  So let's make sure that works
  # too.
  # Flip the order of alleles here to check that aspect too
  sample_attrs <- unlist(list(
    Locus = "A", ExpectedLength1 = 194, ExpectedLength2 = 162))
  sample_data <- testrds("sample_data_guided_flip.rds")
  sample_summary_expected <- testrds("sample_summary_guided_flip.rds")
  sample_summary <- summarize_sample_guided(
    sample_data, sample_attrs, min_locus_reads = 500)
  expect_equal(sample_summary, sample_summary_expected)
})
