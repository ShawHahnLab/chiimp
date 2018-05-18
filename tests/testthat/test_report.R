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
      groups <- c("   162 bp", "   178 bp", "   182 bp", "   194 bp")
      groups <- factor(groups)
      labels <- c("2", "1", "2", "1")
      expect_equal(plot_data$labels, labels)
      expect_equal(plot_data$groups, groups)
    })
  })

# test tabulate_allele_names ----------------------------------------------

  test_that("tabulate_allele_names produces expected data frame", {
    tbl_known <- data.frame(ID = as.character(1:3), stringsAsFactors = FALSE)
    tbl_known$`A_1` <- c("162-c6933c", "162-c6933c", "182-d679e1")
    tbl_known$`A_2` <- c("194-fc013a", "178-d84dc0", "182-d679e1")
    tbl_known$`B_1` <- c("216-c0f11a", "236-321c79", "220-fb9a92")
    tbl_known$`B_2` <- c("252-27c5bf", "240-2a344f", "236-321c79")
    tbl_known$`1_1` <- c("260-9a01fc", "256-c18a06", "276-ea279a")
    tbl_known$`1_2` <- c("280-74dd46", "284-2b3fab", "280-74dd46")
    tbl_known$`2_1` <- c("250-5dacee", "266-2aa675", "238-6cc8ff")
    tbl_known$`2_2` <- c("318-35b7b6", "266-2aa675", "342-2e88c0")
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
