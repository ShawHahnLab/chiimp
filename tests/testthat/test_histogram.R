context("Test histogram function")

with(test_data, {

  # test histogram ----------------------------------------------------------

  test_that("histogram plots histograms of sequence counts versus lengths", {
    # The basic case.
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    sample_data <- results_summary_data$results$samples[["3-2"]]
    output <- histogram(seq_data = seq_data,
                       sample_data = sample_data,
                       cutoff_fraction = 0.05,
                       render = FALSE)
    # The return value of histogram() should be a list of data frames organized
    # by category.
    output_expected <- within(list(), {
      allele <- data.frame(Length = c(318, 350),
                           Count = c(2288, 1270))
      topknown <- data.frame(Length = c(314, 318, 346, 350),
                             Count = c(804, 2288, 138, 1270),
                             SeqName = c("314-ce2338", "318-35b7b6",
                                         "346-b05233", "350-4acdbb"),
                             stringsAsFactors = FALSE)
      topcounts <- data.frame(Length = c(314, 318, 346, 350),
                              Count = c(804, 2288, 138, 1270))
      filt <- data.frame(Length = c(314, 318, 346, 350),
                         Count = c(804, 2288, 138, 1270))
      orig <- data.frame(Length = c(41, 43, 44, 49, 50, 52, 57, 314, 318,
                                    346, 350),
                         Count = c(100, 50, 50, 100, 50, 50, 100, 804, 2288,
                                   138, 1270))
    })
    output_expected <- lapply(output_expected, function(x) {
      x$Length <- as.integer(x$Length)
      x$Count <- as.integer(x$Count)
      x
    })
    expect_equal(output, output_expected)
  })

  test_that("histogram plots histogram with just seq_data provided", {
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    output <- histogram(seq_data = seq_data,
                       render = FALSE)
    output_expected <- within(list(), {
      orig <- data.frame(Length = as.integer(c(41, 43, 44, 49, 50, 52, 57, 314,
                                               318, 346, 350)),
                         Count = as.integer(c(100, 50, 50, 100, 50, 50, 100,
                                              804, 2288, 138, 1270)))
    })
    expect_equal(output, output_expected)
  })

  test_that("histogram works with empty seq_data", {
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]][0, ]
    output <- histogram(seq_data = seq_data,
                       render = FALSE)
    output_expected <- within(list(), {
      orig <- data.frame(Length = as.integer(NA), Count = as.integer(NA))[0, ]
    })
    expect_equal(output, output_expected)
  })

  test_that("histogram plots histogram with empty sample_data", {
    # What if we have sequences present but they all get filtered out for a
    # given sample?  The function should still work.
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    sample_data <- results_summary_data$results$samples[["3-2"]][0, ]
    output <- histogram(seq_data = seq_data,
                       sample_data = sample_data,
                       cutoff_fraction = 0.05,
                       render = FALSE)
    # With sample_data data empty but still provided we should still have the
    # same structure of output, just with zero-row data frames in most of the
    # list entries.
    output_expected <- within(list(), {
      allele <- data.frame(Length = NA, Count = NA)[0, ]
      topknown <- data.frame(Length = NA,
                             Count = NA,
                             SeqName = as.character(NA),
                             stringsAsFactors = FALSE)[0, ]
      topcounts <- data.frame(Length = NA,
                              Count = NA)[0, ]
      filt <- data.frame(Length = NA,
                         Count = NA)[0, ]
      orig <- data.frame(Length = c(41, 43, 44, 49, 50, 52, 57, 314, 318,
                                    346, 350),
                         Count = c(100, 50, 50, 100, 50, 50, 100, 804, 2288,
                                   138, 1270))
    })
    output_expected <- lapply(output_expected, function(x) {
      x$Length <- as.integer(x$Length)
      x$Count <- as.integer(x$Count)
      x
    })
    expect_equal(output, output_expected)
  })

})
