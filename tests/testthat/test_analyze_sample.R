context("Test sample analysis")


test_that("analyze_sample filters and categorizes locus-specific sequences", {
  # analyze_sample should filter a sequence table to entries matching the
  # primer(s), length range, and motif of a specific locus, and assign a
  # category (allele, etc.) to each row.
  seq_data <- readRDS(test_path("data", "analyze_sample", "seq_data.rds"))
  sample_data_expected <- readRDS(test_path(
    "data", "analyze_sample", "sample_data.rds"))
  sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
  expect_equal(sample_data, sample_data_expected)
})
