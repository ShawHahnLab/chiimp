testrds <- function(fname) readRDS(test_path("data", "report", fname))


# test plot_alignment -----------------------------------------------------


test_that("plot_alignment plots alignments for a locus", {
  alignments <- testrds("alignments.rds")
  fp_img <- tempfile()
  png(fp_img)
  plot_data <- plot_alignment(alignments[["A"]])
  dev.off()
  unlink(x = fp_img)
  groups <- c("   162 bp", "   178 bp", "   182 bp", "   194 bp")
  groups <- factor(groups)
  labels <- c("2", "1", "2", "1")
  expect_equal(plot_data$labels, labels)
  expect_equal(plot_data$groups, groups)
})


# test plot_dist_mat ------------------------------------------------------


test_that("plot_dist_mat plots heatmap of distance matrix", {
  dist_mat <- testrds("dist_mat.rds")
  fp_img <- tempfile()
  png(fp_img)
  plot_data <- plot_dist_mat(dist_mat)
  dev.off()
  unlink(x = fp_img)
  # check a few basic things on the pheatmap object
  expect_equal(class(plot_data), "pheatmap")
  expect_equal(
    plot_data$gtable$layout$name,
    c("col_tree", "matrix", "col_names", "row_names", "legend"))
})

test_that("plot_dist_mat can take pheatmap arguments", {
  dist_mat <- testrds("dist_mat.rds")
  fp_img <- tempfile()
  png(fp_img)
  annotations <- data.frame(
    Attribute = LETTERS[seq_along(rownames(dist_mat))])
  rownames(annotations) <- rownames(dist_mat)
  plot_data <- plot_dist_mat(dist_mat, annotation_col = annotations)
  dev.off()
  unlink(x = fp_img)
  # we should still get a pheatmap object back, but this time with the
  # annotation info built into the plot data
  expect_equal(class(plot_data), "pheatmap")
  expect_true("col_annotation" %in% plot_data$gtable$layout$name)
})


# test tabulate_allele_names ----------------------------------------------


test_that("tabulate_allele_names produces expected data frame", {
  results_summary <- testrds("results_summary.rds")
  tbl_known <- testrds("allele_summary.rds")
  tbl <- tabulate_allele_names(results_summary)
  expect_equivalent(tbl, tbl_known)
})

test_that("tabulate_allele_names works on completely blank set", {
  # What if we have a dataset where nothing worked?  The table should be
  # created but should be completely blank.
  tbl_known <- testrds("allele_summary.rds")
  results_summary <- testrds("results_summary.rds")
  for (locus in c("A", "B", "1", "2")) {
    tbl_known[paste(locus, c("1", "2"), sep = "_")] <- ""
  }
  results_summary[, c("Allele1Name", "Allele2Name")] <- NA
  tbl <- tabulate_allele_names(results_summary)
  expect_equivalent(tbl, tbl_known)
})

test_that("tabulate_allele_names keeps selected columns", {
  results_summary <- testrds("results_summary.rds")
  tbl_known <- data.frame(ID = as.character(1:3),
                          Sample = as.character(1:3),
                          Replicate = integer(3) * NA,
                          stringsAsFactors = FALSE)
  # The basic case: sample and replicate columns
  tbl <- tabulate_allele_names(results_summary,
                               extra_cols = c("Sample", "Replicate"))
  expect_true(all(c("Sample", "Replicate") %in% colnames(tbl)))
  expect_equivalent(tbl[, 1:3], tbl_known)
  # What about just one column name?  (We'd better be using drop=FALSE!)
  tbl <- tabulate_allele_names(results_summary,
                               extra_cols = "Sample")
  expect_true("Sample" %in% colnames(tbl))
  expect_false("Replicate" %in% colnames(tbl))
  expect_equivalent(tbl[, 1:2], tbl_known[, 1:2])
  # Unknown column name
  expect_error({
    tabulate_allele_names(results_summary,
                          extra_cols = "DoesNotExist")
  }, "undefined extra columns")
  # Empty column name vector -> same as not given
  tbl <- tabulate_allele_names(results_summary, extra_cols = character())
  expect_true(all(! c("Sample", "Replicate") %in% colnames(tbl)))
})


test_that("tabulate_allele_names sorts rows using extra cols given", {
  # By default tabulate_allele_names sorts in whatever order is given, but
  # there are some edge cases for certain combinations.  The optional
  # extra_cols argument can help enforce a particular output ordering.

  # We'll make the the data have sample 1 replicate 1, sample 2 replicate 1,
  # and then sample 1 replicate 2.
  results_summary <- testrds("results_summary.rds")
  results_summary$Replicate <- "1"
  results_summary$Replicate[results_summary$Sample == "3"] <- "2"
  results_summary$Sample <- rep(c("1", "2", "1"), 4)
  tbl_known <- data.frame(
    ID = c("1-1", "2-1", "1-2"),
    Sample = c("1", "2", "1"),
    Replicate = c("1", "1", "2"),
    stringsAsFactors = FALSE)
  tbl_known$`A_1` <- c("162-c6933c", "162-c6933c", "182-d679e1")
  tbl_known$`A_2` <- c("194-fc013a", "178-d84dc0", "182-d679e1")
  tbl_known$`B_1` <- c("216-c0f11a", "236-321c79", "220-fb9a92")
  tbl_known$`B_2` <- c("252-27c5bf", "240-2a344f", "236-321c79")
  tbl_known$`1_1` <- c("260-9a01fc", "256-c18a06", "276-ea279a")
  tbl_known$`1_2` <- c("280-74dd46", "284-2b3fab", "280-74dd46")
  tbl_known$`2_1` <- c("250-5dacee", "266-2aa675", "238-6cc8ff")
  tbl_known$`2_2` <- c("318-35b7b6", "266-2aa675", "342-2e88c0")

  # Generally, order going in -> order going out, so this leaves sample 1
  # split to either side of sample 2's row.
  tbl <- tabulate_allele_names(results_summary)
  expect_equivalent(tbl, subset(tbl_known, select = -c(Sample, Replicate)))
  # but with extra_cols we can sort while tabulating
  tbl <- tabulate_allele_names(
    results_summary,
    extra_cols = c("Sample", "Replicate"))
  expect_equivalent(tbl, tbl_known[order_entries(tbl_known), ])
  # If we start with sorted input, we can generally expect sorted output
  tbl <- tabulate_allele_names(
    results_summary[order_entries(results_summary), ])
  expect_equivalent(
    tbl,
    subset(
      tbl_known[order_entries(tbl_known), ],
      select = -c(Sample, Replicate)))

  # Partial grid of samples/loci also works, though this is where an edge case
  # comes up about sorting with/without a Locus column included (with, for the
  # input; without, for the output)
  results_summary <- subset(
    results_summary,
    paste(Sample, Replicate, Locus) %in% c(
      "1 1 A", "1 1 B", "2 1 A", "2 1 B", "1 2 B"))
  tbl_known <- tbl_known[
    , c("ID", "Sample", "Replicate", "A_1", "A_2", "B_1", "B_2")]
  tbl_known[3, 4:5] <- NA
  # As before, we get equivalent output sorting as input, generally
  tbl <- tabulate_allele_names(results_summary)
  expect_equivalent(tbl, subset(tbl_known, select = -c(Sample, Replicate)))
  # But, with an incomplete grid of samples/loci, sorting doesn't work like
  # above and we still get the same output.
  tbl <- tabulate_allele_names(
    results_summary[order_entries(results_summary), ])
  expect_equivalent(tbl, subset(tbl_known, select = -c(Sample, Replicate)))
  # Giving extra columns enforces output sorting as expected
  tbl <- tabulate_allele_names(
    results_summary,
    extra_cols = c("Sample", "Replicate"))
  expect_equivalent(tbl, tbl_known[order_entries(tbl_known), ])
  tbl <- tabulate_allele_names(
    results_summary[order_entries(results_summary), ],
    extra_cols = c("Sample", "Replicate"))
  expect_equivalent(tbl, tbl_known[order_entries(tbl_known), ])
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
  results <- list(summary = testrds("results_summary.rds"))
  fp_img <- tempfile()
  png(fp_img)
  plot_data <- plot_heatmap(results, "Stutter")
  dev.off()
  unlink(x = fp_img)
  expect_equal(class(plot_data), "pheatmap")
})

test_that("plot_heatmap handles empty results", {
  # heatmap rendering should still work for a completely empty dataset (all NA
  # entries)
  results <- list(summary = testrds("results_summary.rds"))
  fp_img <- tempfile()
  # empty out certain columns
  results$summary <- zero_summary(results$summary)
  # the function should still work as before
  png(fp_img)
  plot_data <- plot_heatmap(results, "Stutter")
  dev.off()
  unlink(x = fp_img)
  expect_equal(class(plot_data), "pheatmap")
})

test_that("plot_heatmap handles single-value case", {
  # heatmap rendering should still work for a dataset with only one unique
  # value
  results <- list(summary = testrds("results_summary.rds"))
  fp_img <- tempfile()
  # force all entries to a single value
  results$summary$Stutter <- TRUE
  png(fp_img)
  plot_data <- plot_heatmap(results, "Stutter")
  dev.off()
  unlink(x = fp_img)
  expect_equal(class(plot_data), "pheatmap")
})

test_that("plot_heatmap handles single-value case with blanks", {
  results <- list(summary = testrds("results_summary.rds"))
  # heatmap rendering should still work for a dataset with only one unique
  # value plus some NA entries
  fp_img <- tempfile()
  # force all entries to a single value
  results$summary$Stutter <- TRUE
  results$summary[1:4, ] <- zero_summary(results$summary[1:4, ])
  png(fp_img)
  plot_data <- plot_heatmap(results, "Stutter")
  dev.off()
  unlink(x = fp_img)
  expect_equal(class(plot_data), "pheatmap")
})


# test plot_cts_per_locus -------------------------------------------------


test_that("plot_cts_per_locus plots heatmap of counts per matched locus", {
  # It doesn't return anything useful right now, but it should run without
  # errors, at least.
  results <- testrds("results.rds")
  output <- plot_cts_per_locus(results$cts_per_locus, render = FALSE)
  expect_null(output)
})

test_that("plot_cts_per_locus works on completely blank set", {
  # It should still run without errors
  results <- testrds("results.rds")
  results$cts_per_locus[, ] <- 0
  output <- plot_cts_per_locus(results$cts_per_locus, render = FALSE)
  expect_null(output)
})


# test report_genotypes ---------------------------------------------------


test_that("report_genotypes produces expected data frame", {
  # Basic test
  # Largely just a wrapper around tabulate_allele_names, but with a few
  # additional features like NA handling for specific kinds of columns
  tbl_known <- data.frame(
    Sample = as.character(1:3), stringsAsFactors = FALSE)
  tbl_known$`A_1` <- c("162-c6933c", "162-c6933c", "182-d679e1")
  tbl_known$`A_2` <- c("194-fc013a", "178-d84dc0", "182-d679e1")
  tbl_known$`B_1` <- c("216-c0f11a", "236-321c79", "220-fb9a92")
  tbl_known$`B_2` <- c("252-27c5bf", "240-2a344f", "236-321c79")
  tbl_known$`1_1` <- c("260-9a01fc", "256-c18a06", "276-ea279a")
  tbl_known$`1_2` <- c("280-74dd46", "284-2b3fab", "280-74dd46")
  tbl_known$`2_1` <- c("250-5dacee", "266-2aa675", "238-6cc8ff")
  tbl_known$`2_2` <- c("318-35b7b6", "266-2aa675", "342-2e88c0")
  results <- testrds("results.rds")
  tbl <- report_genotypes(results)
  expect_equivalent(tbl, tbl_known)
  # Missing combinations in the input should give blanks in the table, and by
  # default so should missing results
  # (Sample 3 Locus A not given; Sample 2 Locus A blank results)
  results$summary[
    11, c("Allele1Seq", "Allele2Seq", "Allele1Name", "Allele2Name")] <- NA
  results$summary <- results$summary[-12, ]
  tbl <- report_genotypes(results)
  tbl_known[2:3, 8:9] <- ""
  expect_equivalent(tbl, tbl_known)
  # Blanks are the default for NA here but we can specify something else
  tbl <- report_genotypes(results, na.alleles = "X")
  expect_identical(tbl[, 8], c("250-5dacee", "", "X"))
  expect_identical(tbl[, 9], c("318-35b7b6", "", "X"))
})

test_that("report_genotypes handles replicates including NA", {
  # Test for na.replicates argument
  results <- testrds("results.rds")
  # Explicitly label Sample 1 with a replicate, which will make that column
  # show up in the output
  results$summary$Replicate[results$summary$Sample == 1] <- 1
  tbl <- report_genotypes(results)
  expect_identical(tbl$Replicate, c("1", "", ""))
  # Blanks are the default for NA here but we can specify something else
  tbl <- report_genotypes(results, na.replicates = "X")
  expect_identical(tbl$Replicate, c("1", "X", "X"))
})

test_that("report_genotypes uses text for absent sample/locus combos", {
  # Test for na.alleles argument
  # remove one tested combo from the results
  results <- testrds("results.rds")
  results$summary <- subset(results$summary, ! (Sample == 3 & Locus == 2))
  results$files <- results$files[results$summary$Filename]
  results$samples <- results$samples[rownames(results$summary)]
  # by default, an empty string is shown for missing info, indistinguishable
  # from blank results.  Locus 1 should be unaffected, but we should see a
  # blank for sample 3 in Locus 2's first column.
  tbl <- report_genotypes(results)
  expect_equal(tbl[["1_2"]], c("280-74dd46", "284-2b3fab", "280-74dd46"))
  expect_equal(tbl[["2_1"]], c("250-5dacee", "266-2aa675", ""))
  # If we give an na.alleles argument we should be able to get different
  # placeholder text there.
  tbl <- report_genotypes(results, na.alleles = "X")
  expect_equal(tbl[["1_2"]], c("280-74dd46", "284-2b3fab", "280-74dd46"))
  expect_equal(tbl[["2_1"]], c("250-5dacee", "266-2aa675", "X"))
  # That placeholder text should only be applied to allele columns,
  # not elsewhere like Replicate or known ID info columns
  results$summary$Replicate <- rep(1, nrow(results$summary))
  results$summary$Replicate[results$summary$Sample == 3] <- NA
  tbl <- report_genotypes(results)
  expect_equal(tbl$Replicate, c("1", "1", ""))
  tbl <- report_genotypes(results, na.alleles = "X")
  expect_equal(tbl$Replicate, c("1", "1", ""))
  # That's somewhat a special case, though, since Replicate has some
  # NA-handling logic of its own.  How about the identity columns, if present?
  # (Faking the output from find_closest_matches here: nobody has a close
  # match except for sample 3, which matches Bob perfectly)
  closest <- lapply(rownames(tbl), function(entryname) numeric())
  names(closest) <- rownames(tbl)
  closest[["3"]] <- c(Bob = 0)
  tbl <- report_genotypes(results, closest = closest)
  expect_equal(tbl[["Distance"]], c("", "", "0"))
  expect_equal(tbl[["Name"]], c("", "", "Bob"))
  tbl <- report_genotypes(results, closest = closest, na.alleles = "X")
  expect_equal(tbl[["Distance"]], c("", "", "0"))
  expect_equal(tbl[["Name"]], c("", "", "Bob"))
})
