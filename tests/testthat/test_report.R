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

  test_that("tabulate_allele_names works on completely blank set", {
    # What if we have a dataset where nothing worked?  The table should be
    # created but should be completely blank.
    tbl_known <- data.frame(ID = as.character(1:3), stringsAsFactors = FALSE)
    for (locus in c("A", "B", "1", "2")) {
      tbl_known[paste(locus, c("1", "2"), sep = "_")] <- ""
    }
    with(results_summary_data, {
      results$summary[, c("Allele1Name", "Allele2Name")] <- NA
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


# test plot_heatmap -------------------------------------------------------

  # empty out columns of a results summary data frame
  zero_summary <- function(results_summary) {
    for (x in c("Seq", "Count", "Length", "Name")) {
      results_summary[, paste0(c("Allele1", "Allele2"), x)] <- NA
    }
    for (x in c("Homozygous", "Ambiguous", "Stutter", "Artifact")) {
      results_summary[[x]] <- FALSE
    }
    for (x in c("ProminentSeqs")) {
      results_summary[[x]] <- 0
    }
    return(results_summary)
  }

  test_that("plot_heatmap renders heatmap of attribute", {
    # basic test of plot_heatmap.  It should return a pheatmap object.
    with(results_summary_data, {
      fp_img <- tempfile()
      png(fp_img)
      plot_data <- plot_heatmap(results, "Stutter")
      dev.off()
      expect_equal(class(plot_data), "pheatmap")
    })
  })

  test_that("plot_heatmap handles empty results", {
    # heatmap rendering should still work for a completely empty dataset (all NA
    # entries)

    with(results_summary_data, {
      fp_img <- tempfile()
      # empty out certain columns
      results$summary <- zero_summary(results$summary)
      # the function should still work as before
      png(fp_img)
      plot_data <- plot_heatmap(results, "Stutter")
      dev.off()
      expect_equal(class(plot_data), "pheatmap")
    })
  })

  test_that("plot_heatmap handles single-value case", {
    # heatmap rendering should still work for a dataset with only one unique
    # value
    with(results_summary_data, {
      fp_img <- tempfile()
      # force all entries to a single value
      results$summary$Stutter <- TRUE
      png(fp_img)
      plot_data <- plot_heatmap(results, "Stutter")
      dev.off()
      expect_equal(class(plot_data), "pheatmap")
    })
  })

  test_that("plot_heatmap handles single-value case with blanks", {
    # heatmap rendering should still work for a dataset with only one unique
    # value plus some NA entries
    with(results_summary_data, {
      fp_img <- tempfile()
      # force all entries to a single value
      results$summary$Stutter <- TRUE
      results$summary[1:4, ] <- zero_summary(results$summary[1:4, ])
      png(fp_img)
      plot_data <- plot_heatmap(results, "Stutter")
      dev.off()
      expect_equal(class(plot_data), "pheatmap")
    })
  })

# test plot_cts_per_locus -------------------------------------------------


  test_that("plot_cts_per_locus plots heatmap of counts per matched locus", {
    # It doesn't return anything useful right now, but it should run without
    # errors, at least.
    with(results_summary_data, {
      results <- summarize_dataset(results)
      output <- plot_cts_per_locus(results$cts_per_locus, render = FALSE)
      expect_null(output)
    })
  })

  test_that("plot_cts_per_locus works on completely blank set", {
    # It should still run without errors
    with(results_summary_data, {
      results <- summarize_dataset(results)
      results$cts_per_locus[, ] <- 0
      output <- plot_cts_per_locus(results$cts_per_locus, render = FALSE)
      expect_null(output)
    })
  })

})
