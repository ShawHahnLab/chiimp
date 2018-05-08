context("Test reporting functions")

with(test_data, {

# test plot_alignment -----------------------------------------------------

  test_that("plot_alignment plots alignments for a locus", {
    with(results_summary_data, {
      alignments <- align_alleles(results$summary)
      fp_img <- tempfile()
      png(fp_img)
      plot_data <- plot_alignment(alignments[["A"]])
      dev.off()
      groups <- c("   162 bp",
                  "   174 bp",
                  "   178 bp",
                  "   182 bp",
                  "   194 bp")
      groups <- factor(groups)
      labels <- c("1", "2", "1", "1", "1")
      expect_equal(plot_data$labels, labels)
      expect_equal(plot_data$groups, groups)
    })
  })



# test histogram ----------------------------------------------------------

  test_that("histogram plots histograms of sequence counts versus lengths", {
    # The basic case.
    data_summary <- results_summary_data$results$summary["3-2", ]
    seq_data <- results_summary_data$results$files[[data_summary$Filename]]
    sample_data <- results_summary_data$results$samples[["3-2"]]
    output <- str_hist(seq_data = seq_data,
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
    output <- str_hist(seq_data = seq_data,
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
    output <- str_hist(seq_data = seq_data,
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
    output <- str_hist(seq_data = seq_data,
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

# test tabulate_allele_names ----------------------------------------------

  test_that("tabulate_allele_names produces expected data frame", {
    tbl_known <- data.frame(ID = as.character(1:3), stringsAsFactors = FALSE)
    tbl_known$`A_1` <- c("162-c6933c", "178-d84dc0", "174-8a43ea")
    tbl_known$`A_2` <- c("194-fc013a", "182-d679e1", "174-8a43ea")
    tbl_known$`B_1` <- c("220-fb9a92", "212-6d4afb", "224-053ec2")
    tbl_known$`B_2` <- c("244-3c2ff2", "220-fb9a92", "224-053ec2")
    tbl_known$`1_1` <- c("260-9a01fc", "276-ea279a", "252-a5eee8")
    tbl_known$`1_2` <- c("272-292a2a", "288-201179", "276-ea279a")
    tbl_known$`2_1` <- c("278-ae70f3", "314-ce2338", "318-35b7b6")
    tbl_known$`2_2` <- c("342-2e88c0", "346-b05233", "350-4acdbb")
    with(results_summary_data, {
      tbl <- tabulate_allele_names(results$summary)
      expect_equivalent(tbl, tbl_known)
    })
  })

  test_that("tabulate_allele_names keeps selected columns", {
    tbl_known <- data.frame(ID = as.character(1:3),
                            Sample = as.character(1:3),
                            Replicate = integer(3) * NA,
                            stringsAsFactors = FALSE)
    with(results_summary_data, {
      tbl <- tabulate_allele_names(results$summary,
                                   extra_cols = c("Sample", "Replicate"))
      tbl <- tbl[, 1:3]
      expect_equivalent(tbl, tbl_known)
    })
  })

})
