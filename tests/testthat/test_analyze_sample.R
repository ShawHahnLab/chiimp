context("Test sample analysis")

with(test_data, {

  test_that("analyze_sample categorizes sequences", {
    skip("test not yet implemented")
  })

  test_that("analyze_sample records analysis attributes", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    expect_equal(attr(sample_data, "fraction.min"), 0.05)
  })

})
