context("Test sample summarization")

with(test_data, {

# util --------------------------------------------------------------------

  check.seqs1A_summary <- function(data,
                                   count.locus=4500,
                                   allele1.count=2803,
                                   allele2.count=1300) {
    expect_equal(names(data), sample.summary.cols)
    with(data, {
      expect_equal(Allele1Seq,
                   gsub("[\n ]*", "",
                       "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                        ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                        TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACACAG
                        TTGTGTGAGCCAGTC"))
      expect_equal(Allele2Seq,
                   gsub("[\n ]*", "",
                       "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                        ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                        TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT
                        AGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC"))
      expect_equal(Allele1Count, allele1.count)
      expect_equal(Allele1Length, 162)
      expect_equal(Allele2Count, allele2.count)
      expect_equal(Allele2Length, 194)
      expect_equal(Homozygous, FALSE)
      expect_equal(Stutter, FALSE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, count.locus)
      expect_equal(ProminentSeqs, 3)
    })
  }

# test summarize_sample ---------------------------------------------------

  test_that("summarize_sample summarizes sample attributes", {
    sample.data <- analyze_sample(seqs1$A, locus_attrs, 3)
    sample.summary <- summarize_sample(sample.data, list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    check.seqs1A_summary(sample.summary)
  })

  test_that("summarize_sample handles completely emtpy input sample data", {
    sample.data <- analyze_sample(c(), locus_attrs, 3)
    sample.summary <- summarize_sample(sample.data,  list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    expect_equal(names(sample.summary), sample.summary.cols)
    with(sample.summary, {
      expect_equal(Allele1Seq, as.character(NA))
      expect_equal(Allele1Count, as.integer(NA))
      expect_equal(Allele1Length, as.integer(NA))
      expect_equal(Allele2Seq, as.character(NA))
      expect_equal(Allele2Count, as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous, FALSE)
      expect_equal(Stutter, FALSE)
      expect_equal(CountTotal, 0)
      expect_equal(CountLocus, 0)
      expect_equal(ProminentSeqs, 0)
    })
  })

  test_that("summarize_sample handles empty sequences in input sample data", {
    seqs <- seqs1$A
    seqs[1:100] <- "" # empty out a segment of the vector
    sample.data <- analyze_sample(seqs, locus_attrs, 3)
    sample.summary <- summarize_sample(sample.data,  list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    # Nothing should change in the output, except that we zeroed out 90 reads that
    # would otherwise get counted (the rest were already set to off-target
    # junk during setup).
    check.seqs1A_summary(sample.summary,
                         count.locus = 4414,
                         allele1.count = 2748,
                         allele2.count = 1276)
  })

  test_that("summarize_sample marks stutter removal", {
    sample.data <- analyze_sample(seqs3$A, locus_attrs, 3)
    sample.summary <- summarize_sample(sample.data,  list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    with(sample.summary, {
      expect_equal(Allele1Seq, paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGA",
                                      "TAGATAGATAGATAGATAGATAGATAGATA",
                                      "GATAGATAGATAGATAGATAGATAGATAGA",
                                      "TAGATAGATAGATAGATAGATAGATAGATA",
                                      "GATAGATAGATAGATAGATAGATAGATAGA",
                                      "TAGACACAGTTGTGTGAGCCAGTC"))
      expect_equal(Allele1Count, 3971)
      expect_equal(Allele1Length, 174)
      expect_equal(Allele2Seq, as.character(NA))
      expect_equal(Allele2Count, as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous, TRUE)
      expect_equal(Stutter, TRUE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, 4500)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample handles multiple stutter sequences", {
    # If multiple candidate allele sequences are marked as potential stutter, they
    # should all be skipped, not just the first.
    sample.data <- analyze_sample(seqs3$A, locus_attrs, 3)
    # Replace the third entry with a different stutter sequence.  Munge the counts
    # around to still total correctly.
    tot <- sum(sample.data$Count)
    sample.data[3, ] <- sample.data[2, ]
    sample.data[3, "Seq"] <- sub("TAGA", "TACA", sample.data[3, "Seq"])
    sample.data[3, "Count"] <- 410
    sample.data[3, "FractionOfTotal"] <- 410/tot
    sample.data[3, "FractionOfLocus"] <- 410/tot
    sample.data[4:12, "Count"] <- 10
    sample.data[4:12, "FractionOfTotal"] <- 10/tot
    sample.data[4:12, "FractionOfLocus"] <- 10/tot
    sample.summary <- summarize_sample(sample.data,  list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    with(sample.summary, {
      expect_equal(Allele1Seq, paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGA",
                                      "TAGATAGATAGATAGATAGATAGATAGATA",
                                      "GATAGATAGATAGATAGATAGATAGATAGA",
                                      "TAGATAGATAGATAGATAGATAGATAGATA",
                                      "GATAGATAGATAGATAGATAGATAGATAGA",
                                      "TAGACACAGTTGTGTGAGCCAGTC"))
      expect_equal(Allele1Count, 3971)
      expect_equal(Allele1Length, 174)
      expect_equal(Allele2Seq, as.character(NA))
      expect_equal(Allele2Count, as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous, TRUE)
      expect_equal(Stutter, TRUE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, 4910)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample counts prominent sequences", {
    sample.data.1B       <- analyze_sample(seqs1$B, locus_attrs, 3)
    sample.summary.1B    <- summarize_sample(sample.data.1B,  list(Locus="B"),
                                           fraction.min = 0.05, counts.min = 500)
    sample.data.2B       <- analyze_sample(seqs2$B, locus_attrs, 3)
    sample.summary.2B    <- summarize_sample(sample.data.2B, list(Locus="B"),
                                           fraction.min = 0.05, counts.min = 500)
    sample.data.3B       <- analyze_sample(seqs3$B, locus_attrs, 3)
    sample.summary.3B    <- summarize_sample(sample.data.3B, list(Locus="B"),
                                           fraction.min = 0.05, counts.min = 500)
    sample.data.empty    <- analyze_sample(c(), locus_attrs, 3)
    sample.summary.empty <- summarize_sample(sample.data.empty, list(Locus="B"),
                                           fraction.min = 0.05, counts.min = 500)

    expect_equal(sample.summary.1B$ProminentSeqs,    3)
    expect_equal(sample.summary.2B$ProminentSeqs,    4)
    expect_equal(sample.summary.3B$ProminentSeqs,    1)
    expect_equal(sample.summary.empty$ProminentSeqs, 0)
    # Despite having stutter-y peaks the first two did not have a potential allele
    # removed, so Stutter == FALSE.  The third had stutter large enough to look
    # like an allele so it was removed.  The last one was empty so there was
    # nothing to have stutter from.
    expect_equal(sample.summary.1B$Stutter,    FALSE)
    expect_equal(sample.summary.2B$Stutter,    FALSE)
    expect_equal(sample.summary.3B$Stutter,    TRUE)
    expect_equal(sample.summary.empty$Stutter, FALSE)
  })

  test_that("summarize_sample rejects low-count samples", {
    sample.data <- analyze_sample(seqs1$A, locus_attrs, 3)
    # Here we check that the filtered-counts-thresholding is applied, by forcing
    # the counts to a low number.  This should still report some stats but should
    # leave out the allele1/allele2 information.
    sample.data$Count <- sample.data$Count / 100
    sample.summary <- summarize_sample(sample.data, list(Locus="A"),
                                       fraction.min = 0.05, counts.min = 500)
    with(sample.summary, {
      expect_equal(Allele1Seq, as.character(NA))
      expect_equal(Allele1Count, as.integer(NA))
      expect_equal(Allele1Length, as.integer(NA))
      expect_equal(Allele2Seq, as.character(NA))
      expect_equal(Allele2Count, as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous, FALSE)
      expect_equal(Stutter, FALSE)
      expect_equal(CountTotal, 50)
      expect_equal(CountLocus, 45)
      expect_equal(ProminentSeqs, 3)
    })
  })

  test_that("summarize_sample warns of missing locus name", {
    # summarize_sample() should be able to tell if an invalid (as per the earlier
    # processing) locus name is given, because it won't be in the levels of the
    # MatchingLocus factor of the data frame from analyze_sample().  I think this
    # would only come about when locus names given by prepare.dataset() don't
    # match what's in locus_attrs.
    skip("test not yet implemented")
  })

})
