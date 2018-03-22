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


# test tabulate_allele_names ----------------------------------------------

  test_that("tabulate_allele_names produces expected data frame", {
    fail("not implemented")
    with(results_summary_data, {
      tbl <- tabulate_allele_names(results$summary)
      # TODO clean up stuff betwen analyze_dataset and summarize_dataset in
      # chiimp.R and put the factoring/ordering code inside analyze_dataset (?)
    })
  })

})
