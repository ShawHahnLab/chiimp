testrds <- function(fname) readRDS(test_path("data", "analyze_seqs", fname))


# test analyze_seqs -------------------------------------------------------


test_that("analyze_seqs tabulates sequences", {
  seqs <- testrds("seqs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data.rds")
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs handles completely emtpy input vector", {
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_empty.rds")
  seq_data <- analyze_seqs(c(), locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs handles empty sequences", {
  # As in seqs.rds but with the first hundred sequences set to empty strings
  seqs <- testrds("seqs_stubs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_stubs.rds")
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs can use reverse primers", {
  seqs <- testrds("seqs.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_rev_primers.rds")
  # As is, everything still works if we enable the use of reverse primers from
  # the locus attributes table.
  seq_data <- analyze_seqs(seqs, locus_attrs, 3, use_reverse_primers = TRUE)
  expect_equal(seq_data, seq_data_expected)
  # Now, try with a reverse primer that won't match up.  Here A's reverse
  # primer is replaced with B's.
  locus_attrs_mod <- testrds("locus_attrs_rev_primer_mod.rds")
  seq_data_expected_mod <- testrds("seq_data_rev_primers_mod.rds")
  # No more locus A matched since we replaced the reverse primer with B's in
  # the locus attributes but the forward primer observed still matches A's.
  seq_data <- analyze_seqs(seqs, locus_attrs_mod, 3, use_reverse_primers = TRUE)
  expect_equal(seq_data, seq_data_expected_mod)
})

test_that("analyze_seqs can use reverse primers and auto-revcmp", {
  seqs <- testrds("seqs.rds")
  # If we supply the reverse primers in their orientation on R2 (reverse
  # complement of what's there in the other examples), it should still work as
  # expected so long as we specify reverse_primer_r1 = FALSE.
  locus_attrs_mod <- testrds("locus_attrs_rev_primer_revcmp.rds")
  seq_data_expected <- testrds("seq_data_rev_primers.rds")
  seq_data <- analyze_seqs(
    seqs, locus_attrs_mod, 3,
    use_reverse_primers = TRUE,
    reverse_primer_r1 = FALSE)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs does not mark high-count entries as stutter", {
  # In this version I've  replaced a few high-count entries with one that was
  # considered stutter in the above test.  Now none of those top sequences
  # should be considered stutter.
  seqs <- testrds("seqs_stutter_filter_check.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_stutter_filter_check.rds")
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs works with varied threshold for stutter counts", {
  # Here I've removed seqs at 158 nt and 54 nt and put those counts into the
  # 190 nt seq instead, making that more abundant compared to 194 (ratio of
  # 0.34).  With the default ratio setting of 1/3 it would no longer be
  # considered stutter.
  seqs <- testrds("seqs_stutter_threshold_check.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_stutter_threshold_orig.rds")
  seq_data_expected_mod <- testrds("seq_data_stutter_threshold_mod.rds")
  # With default settings the 190 nt seq is no longer considered stutter.
  # With a higher ratio_max value it still is
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
  seq_data_mod <- analyze_seqs(
    seqs, locus_attrs, 3, stutter.count.ratio_max = 1 / 2)
  expect_equal(seq_data_mod, seq_data_expected_mod)
})

test_that("analyze_seqs marks artifact entries", {
  # Here I took the first stutter entry in the other cases and made it a
  # different-by-one-nt version of the most abundant sequence
  seqs <- testrds("seqs_artifact.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_artifact.rds")
  # The third entry should be marked an artifact of the first
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})

test_that("analyze_seqs marks ambiguous entries", {
  # sequences that contain non-ACTG characters should be marked TRUE in the
  # Ambiguous column (interpreting those as "N" or similar).
  # Here I substituted an N at a specific position for one sequence's entries
  seqs <- testrds("seqs_ambig.rds")
  locus_attrs <- testrds("locus_attrs.rds")
  seq_data_expected <- testrds("seq_data_ambig.rds")
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  expect_equal(seq_data, seq_data_expected)
})
