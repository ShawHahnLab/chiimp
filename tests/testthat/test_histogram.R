testrds <- function(fname) readRDS(test_path("data", "histogram", fname))


# test histogram ----------------------------------------------------------


test_that("histogram plots histograms of sequence counts versus lengths", {
  seq_data <- testrds("seq_data.rds")
  sample_data <- testrds("sample_data.rds")
  png(fp_devnull)
  output <- histogram(
    seq_data = seq_data, sample_data = sample_data, cutoff_fraction = 0.05)
  dev.off()
  # The return value of histogram() should be a list of data frames organized
  # by category.
  output_expected <- testrds("output.rds")
  expect_equal(output, output_expected)
})

test_that("histogram plots histogram with just seq_data provided", {
  seq_data <- testrds("seq_data.rds")
  png(fp_devnull)
  output <- histogram(seq_data = seq_data)
  dev.off()
  output_expected <- testrds("output_seq_data_only.rds")
  expect_equal(output, output_expected)
})

test_that("histogram works with empty seq_data", {
  seq_data <- testrds("seq_data.rds")[0, ]
  png(fp_devnull)
  expect_silent(
    output <- histogram(seq_data = seq_data))
  dev.off()
  output_expected <- within(list(), {
    orig <- data.frame(Length = as.integer(NA), Count = as.integer(NA))[0, ]
  })
  expect_equal(output, output_expected)
})

test_that("histogram plots histogram with empty sample_data", {
  # What if we have sequences present but they all get filtered out for a
  # given sample?  The function should still work.
  seq_data <- testrds("seq_data.rds")
  sample_data <- testrds("sample_data.rds")[0, ]
  png(fp_devnull)
  expect_silent(
    output <- histogram(
      seq_data = seq_data, sample_data = sample_data, cutoff_fraction = 0.05))
  dev.off()
  output_expected <- testrds("output_empty_sample_data.rds")
  expect_equal(output, output_expected)
})
