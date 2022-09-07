context("Test histogram function")

with(test_data, {

  # test histogram ----------------------------------------------------------

  test_that("histogram plots histograms of sequence counts versus lengths", {
    # The basic case.
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    sample_data <- results_summary_data$results$samples[["3-2"]]
    png(fp_devnull)
    output <- histogram(seq_data = seq_data,
                       sample_data = sample_data,
                       cutoff_fraction = 0.05)
    dev.off()
    # The return value of histogram() should be a list of data frames organized
    # by category.
    output_expected <- within(list(), {
      allele <- data.frame(Length = c(342, 238),
                           Count = c(2615, 1217))
      topknown <- data.frame(Length = c(238, 342),
                             Count = c(1217, 2615),
                             SeqName = c("238-6cc8ff", "342-2e88c0"),
                             stringsAsFactors = FALSE)
      topcounts <- data.frame(Length = c(234, 238, 338, 342),
                              Count = c(99, 1217, 535, 2615))
      filt <- data.frame(Length = c(234, 238, 338, 342),
                         Count = c(99, 1217, 535, 2615))
      orig <- data.frame(Length = c(40, 41, 43, 46, 47, 48, 50, 162, 166, 220,
                                    224, 232, 234, 236, 238, 264, 268, 272, 338,
                                    342),
                         Count = c(100, 50, 50, 150, 50, 33, 50, 4, 13, 2, 10,
                                   1, 99, 4, 1217, 2, 11, 4, 535, 2615))
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
    png(fp_devnull)
    output <- histogram(seq_data = seq_data)
    dev.off()
    output_expected <- within(list(), {
      orig <- data.frame(Length = as.integer(c(40, 41, 43, 46, 47, 48, 50, 162,
                                               166, 220, 224, 232, 234, 236,
                                               238, 264, 268, 272, 338, 342)),
                         Count = as.integer(c(100, 50, 50, 150, 50, 33, 50, 4,
                                              13, 2, 10, 1, 99, 4, 1217, 2, 11,
                                              4, 535, 2615)))
    })
    expect_equal(output, output_expected)
  })

  test_that("histogram works with empty seq_data", {
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]][0, ]
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
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    sample_data <- results_summary_data$results$samples[["3-2"]][0, ]
    png(fp_devnull)
    expect_silent(
      output <- histogram(seq_data = seq_data,
                          sample_data = sample_data,
                          cutoff_fraction = 0.05))
    dev.off()
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
      orig <- data.frame(Length = as.integer(c(40, 41, 43, 46, 47, 48, 50, 162,
                                               166, 220, 224, 232, 234, 236,
                                               238, 264, 268, 272, 338, 342)),
                         Count = as.integer(c(100, 50, 50, 150, 50, 33, 50, 4,
                                              13, 2, 10, 1, 99, 4, 1217, 2, 11,
                                              4, 535, 2615)))
    })
    output_expected <- lapply(output_expected, function(x) {
      x$Length <- as.integer(x$Length)
      x$Count <- as.integer(x$Count)
      x
    })
    expect_equal(output, output_expected)
  })

})
