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
      groups <- c("   162 bp", "   174 bp", "   178 bp", "   182 bp", "   194 bp")
      groups <- factor(groups)
      labels <- c("1", "2", "1", "1", "1")
      expect_equal(plot_data$labels, labels)
      expect_equal(plot_data$groups, groups)
    })
  })



# test histogram ----------------------------------------------------------

  test_that("histogram plots histograms of sequence counts versus lengths", {

    data <- results_summary_data$results$data[["3-2"]]
    data_summary <- results_summary_data$results$summary["3-2", ]

    # The return value of histogram() should be a list of data frames organized
    # by category.  Attributes will define colors and limits and such on both
    # the data frames and the whole object.
    # TODO finish defining these.  It should be all info needed to execute the
    # actual plot commands.  (Am I just re-creating the concept of ggplot?)
    output_expected <- within(list(), {
      total <- data.frame(Count  = c(2288, 1270, 804, 138, 100, 50,
                                     50,     50,  50,  50,  50, 50, 50),
                          Length = c(318, 350, 314, 346, 41, 50,
                                      44,  57,  52,  49, 43, 57, 49))
      filtered <- data.frame(Count  = c(2288, 1270, 804, 138),
                             Length = c(318, 350, 314, 346))
      unique <- data.frame(Count  = c(2288, 1270, 804, 138),
                           Length = c(318, 350, 314, 346))
      known <- data.frame(Count   = c(2288, 1270, 804, 138),
                          Length  = c(318, 350, 314, 346),
                          SeqName = c("318-35b7b6", "350-4acdbb",
                                      "314-ce2338", "346-b05233"))
      alleles <- data.frame(Count  = c(2288, 1270),
                            Length = c(318, 350))
    })

    attr(output_expected$total,    "color") <- "#000000FF"
    attr(output_expected$filtered, "color") <- "#888888FF"
    attr(output_expected$unique,   "color") <- "#FFAAAAFF"
    attr(output_expected$known,    "color") <- "#0000FFFF"
    attr(output_expected$alleles,  "color") <- "#FF0000FF"
    cutoff <- 250
    xlim.filt <- c(314, 350)
    ymax <- 2288
    attr(output_expected, "cutoff") <- list(h = cutoff,
                                            col = "#00000080")
    attr(output_expected, "region") <- list(x = rep(xlim.filt, each = 2),
                                            y = c(0, ymax, ymax, 0),
                                            col = "#0000001E")

    output <- histogram(samp = data, sample.summary = data_summary,
                        cutoff_fraction = 0.05)

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
