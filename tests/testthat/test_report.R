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
