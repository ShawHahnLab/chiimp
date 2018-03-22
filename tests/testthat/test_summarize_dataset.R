context("Test dataset summarization")

with(test_data, {

# test_summarize_dataset --------------------------------------------------

  test_that("summarize_dataset produces additional summaries", {
    # Just a general check of the overall process and binding these results
    # together in a list; more rigorous checks are below since those are public
    # functions too.
    with(results_summary_data, {
      results_mod <- summarize_dataset(results)
      names_observed <- names(results_mod)
      names_expected <- c(names(results), "alignments", "dist_mat")
      expect_equal(names_observed,
                   names_expected)
    })
  })

# test_make_dist_mat ------------------------------------------------------

  test_that("make_dist_mat produces a valid distance matrix", {
    with(results_summary_data, {
      dist_mat <- make_dist_mat(results$summary)
      # row and column names should both be equal to the sample identifiers
      samps <- unique(dataset$Sample)
      expect_equal(rownames(dist_mat), samps)
      expect_equal(colnames(dist_mat), samps)
      # For our test data samples are far apart from each other, but with no
      # missing entries the diagonal is zero.
      dists <- matrix(7, nrow = length(samps), ncol = length(samps))
      dists[1, 3] <- 8
      dists[3, 1] <- 8
      diag(dists) <- 0
      expect_true(all(dist_mat == dists))
    })
  })

  test_that("make_dist_mat_known produces sample-to-individual distance matrix", {
    skip("test not yet implemented")
  })

# test_calc_genotype_distance ---------------------------------------------

  test_that("calc_genotype_distance scores genotypes correctly", {
    g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
    g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
    d <- calc_genotype_distance(g1, g2)
    expect_equal(d, 4)
  })

  test_that("calc_genotype_distance handles multiple data types", {
    # As written it should work with character, numeric, whatever-- it just uses
    # the != operator.
    g.num <- c(1, 2, 3, 4)
    g.alpha <- c("a", "b", "c", "d")
    d <- calc_genotype_distance(g.num, g.alpha)
    d.num <- calc_genotype_distance(g.num, g.num)
    d.alpha <- calc_genotype_distance(g.alpha, g.alpha)
    expect_equal(d, 4)
    expect_equal(d.num, 0)
    expect_equal(d.alpha, 0)
  })

  test_that("calc_genotype_distance is commutative", {
    g1 <- c(5, 10, 200, 204, 37, 37, 180, 184, 190, 190)
    g2 <- c(5,  5, 204, 208, 37, 39, 176, 180, 190, 190)
    d.fwd <- calc_genotype_distance(g1, g2)
    d.rev <- calc_genotype_distance(g2, g1)
    expect_equal(d.fwd, 4)
    expect_equal(d.rev, 4)
  })

  test_that("calc_genotype_distance handles NA entries", {
    g1 <- c(5,   5, 200, 204, 37, 39)
    g2 <- c(NA, NA, 204, 208, 37, 37)
    g3 <- c(NA, NA, 200, 204, 37, 39)

    # The two NA entries will by default count towards the distance whether the
    # other genotype has an NA there or not.  If na.reject is FALSE, NAs are
    # treated like any other entry.
    d12 <- calc_genotype_distance(g1, g2)
    d12.na <- calc_genotype_distance(g1, g2, na.reject = FALSE)
    d23 <- calc_genotype_distance(g2, g3)
    d23.na <- calc_genotype_distance(g2, g3, na.reject = FALSE)
    expect_equal(d12,    4)
    expect_equal(d12.na, 2)
    expect_equal(d23,    4)
    expect_equal(d23.na, 2)
  })

  test_that("calc_genotype_distance invalid genotype lengths", {
    # These are valid but not the same length
    g1 <- c(5,   5, 200, 204)
    g2 <- c(NA, NA, 204, 208, 37, 39)
    # This has an odd length
    g3 <- c(NA, NA, 204, 208, 37)

    expect_warning(calc_genotype_distance(g1, g2))
    expect_warning(calc_genotype_distance(g3, g3))
  })


# test_align_alleles ------------------------------------------------------

  test_that("align_alleles produces sequence alignments", {
    with(results_summary_data, {
      alignments <- align_alleles(results$summary)
      expect_equal(names(alignments), levels(results$summary$Locus))
      expect_equal(names(as.character(alignments[["A"]])),
                   c("182_1", "194_1", "178_1", "174_2", "162_1"))
      expect_equal(as.character(alignments[["A"]])[[1]],
                   paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGA",
                          "TAGATAGATAGATAGATAGATAGATAGATA",
                          "GATAGATAGATAGATAGATAGATAGATAGA",
                          "TAGATAGATAGATAGATAGATAGATAGATA",
                          "GATAGATAGATAGATAGATAGATAGATAGA",
                          "TAGATAG------------ATAGACACAGT",
                          "TGTGTGAGCCAGTC"))
    })
  })

  test_that("align_alleles produces per-allele alignments", {
    # By default the allele sequences are dereplicated and then aligned, but this
    # tests the other option.
    with(results_summary_data, {
      alignments <- align_alleles(results$summary, derep = FALSE)
      expect_equal(names(alignments), levels(results$summary$Locus))
      n1 <- paste(rownames(subset(results$summary, Locus == "A")), 1, sep = "_")
      n2 <- paste(rownames(subset(results$summary, Locus == "A")), 2, sep = "_")
      expect_equal(sort(names(as.character(alignments[["A"]]))),
                   sort(c(n1, n2)))
    })
  })

  test_that("align_alleles works for empty sequences", {
    # Empty sequences should be handled automatically before the msa function is
    # called.
    with(results_summary_data, {
      # Empty out one entry for locus A
      idx <- which(results$summary$Locus == "A")
      results$summary[idx[1], "Allele1Seq"] <- NA
      results$summary[idx[1], "Allele2Seq"] <- NA
      # Locus A's alignment should be OK, just with one stub entry
      alignments <- align_alleles(results$summary)
      expect_true(any(grepl("^-+$", as.character(alignments$A))))
    })
  })

  test_that("align_alleles works for empty sequence sets", {
    # Completely empty sequence lists for a given locus are a special case; for
    # those we should just get "NA" for the whole entry.
    with(results_summary_data, {
      # Empty out all sequences for locus A
      idx <- results$summary$Locus == "A"
      results$summary[idx, "Allele1Seq"] <- NA
      results$summary[idx, "Allele2Seq"] <- NA
      # This should still work, just with NA for A's alignment
      alignments <- align_alleles(results$summary)
      expect_equal(alignments$A, NULL)
    })
  })

  test_that("align_alleles works for identical sequence sets", {
    # Apparently msa() fails when given a single sequence (derep=TRUE for
    # align_alleles), throwing a "There is an invalid aln file!" error.
    # I'll do its validating for it and account for this inside align_alleles.
    with(results_summary_data, {
      # Overwrite Locus A seqs with same content
      idx <- which(results$summary$Locus == "A")
      results$summary[idx, "Allele1Seq"] <- results$summary[idx[1], "Allele1Seq"]
      results$summary[idx, "Allele2Seq"] <- NA
      alignments <- align_alleles(results$summary)
      # This should still work, just with the same sequence back for A's
      # alignment
      expect_equal(unname(alignments$A),
                   as.character(results$summary[idx[1], "Allele1Seq"]))
    })
  })


# test_tally_cts_per_locus ------------------------------------------------

  test_that("tally_cts_per_locus counts matches per locus per sample", {
    # Right now the function is called in full_analysis, but it really should be
    # rolled into summarize_dataset.
    skip("test not yet implemented")
  })

})
