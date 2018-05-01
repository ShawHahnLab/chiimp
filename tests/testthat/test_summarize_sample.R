context("Test sample summarization")

with(test_data, {

# util --------------------------------------------------------------------

  check.seqs1A_summary <- function(data,
                                   count.locus=4500,
                                   allele1.count=2803,
                                   allele2.count=1300,
                                   ord = 1:2) {
    expect_equal(names(data), sample.summary.cols)
    with(data, {
      alleles <- c(Allele1Seq, Allele2Seq)[ord]
      counts  <- c(Allele1Count, Allele2Count)[ord]
      lengths <- c(Allele1Length, Allele2Length)[ord]
      expect_equal(alleles[1],
                   gsub("[\n ]*", "",
                       "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                        ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                        TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACACAG
                        TTGTGTGAGCCAGTC"))
      expect_equal(alleles[2],
                   gsub("[\n ]*", "",
                       "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                        ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                        TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT
                        AGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC"))
      expect_equal(counts[1],     allele1.count)
      expect_equal(lengths[1],    162)
      expect_equal(counts[2],     allele2.count)
      expect_equal(lengths[2],    194)
      expect_equal(Homozygous,    FALSE)
      expect_equal(Ambiguous,     FALSE)
      expect_equal(Stutter,       FALSE)
      expect_equal(Artifact,      FALSE)
      expect_equal(CountTotal,    5000)
      expect_equal(CountLocus,    count.locus)
      expect_equal(ProminentSeqs, 2)
    })
  }

  expect_empty <- function(sample_summary) {
    with(sample_summary, {
      expect_equal(Allele1Seq,    as.character(NA))
      expect_equal(Allele1Count,  as.integer(NA))
      expect_equal(Allele1Length, as.integer(NA))
      expect_equal(Allele2Seq,    as.character(NA))
      expect_equal(Allele2Count,  as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous,    FALSE)
      expect_equal(Ambiguous,     FALSE)
      expect_equal(Stutter,       FALSE)
      expect_equal(Artifact,      FALSE)
      expect_equal(CountTotal,    0)
      expect_equal(CountLocus,    0)
      expect_equal(ProminentSeqs, 0)
    })
  }

# test summarize_sample ---------------------------------------------------

  test_that("summarize_sample summarizes sample attributes", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"),
                                  fraction.min = 0.05)
    sample_summary <- summarize_sample(sample_data, list(Locus = "A"),
                                       counts.min = 500)
    check.seqs1A_summary(sample_summary)
  })

  test_that("summarize_sample handles completely empty sample data", {
    seq_data <- analyze_seqs(c(), locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data,  list(Locus = "A"),
                                       counts.min = 500)
    expect_equal(names(sample_summary), sample.summary.cols)
    expect_empty(sample_summary)
  })

  test_that("summarize_sample handles empty sequences in input sample data", {
    seqs <- seqs1$A
    seqs[1:100] <- "" # empty out a segment of the vector
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data,  list(Locus = "A"),
                                       counts.min = 500)
    # Nothing should change in the output, except that we zeroed out 90 reads
    # that would otherwise get counted (the rest were already set to off-target
    # junk during setup).
    check.seqs1A_summary(sample_summary,
                         count.locus = 4414,
                         allele1.count = 2748,
                         allele2.count = 1276)
  })

  test_that("summarize_sample marks stutter removal", {
    seq_data <- analyze_seqs(seqs3$A, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data,  list(Locus = "A"),
                                       counts.min = 500)
    with(sample_summary, {
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
      expect_equal(Ambiguous, FALSE)
      expect_equal(Stutter, TRUE)
      expect_equal(Artifact, FALSE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, 4500)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample marks removal of multiple artifact types", {
    # There are separate columns in the output for ambiguous/stutter/artifact
    # sequences and each one can independently be TRUE if a sequence was skipped
    # due to that particular aspect.  (Previously only the first would be
    # marked, and it was incorrectly assumed that a secondary sequence was
    # always the one that may have been filtered.)
    # First take an example sample and split half of the "correct" reads with an
    # ambiguous sequence.  This also leaves some existing stutter reads in
    # place.
    s <- seqs3$A
    idx <- (which(s == s[1]))[c(TRUE, FALSE)] # every other matching index
    s[idx] <- gsub(".$", "N", s[idx]) # replace last character with N
    seq_data <- analyze_seqs(s, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data,  list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data, counts.min = 500)
    with(sample_summary, {
      expect_equal(Homozygous, TRUE)
      expect_equal(Ambiguous, TRUE)
      expect_equal(Stutter, TRUE)
      expect_equal(Artifact, FALSE)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample handles multiple stutter sequences", {
    # If multiple candidate allele sequences are marked as potential stutter,
    # they should all be skipped, not just the first.
    seq_data <- analyze_seqs(seqs3$A, locus_attrs, 3)
    # Replace the third entry with a different stutter sequence.  Munge the
    # counts around to still total correctly.
    tot <- sum(seq_data$Count)
    seq_data[3, ] <- seq_data[2, ]
    seq_data[3, "Seq"] <- sub("TAGA", "TACA", seq_data[3, "Seq"])
    seq_data[3, "Count"] <- 410
    seq_data[3, "FractionOfTotal"] <- 410 / tot
    seq_data[3, "FractionOfLocus"] <- 410 / tot
    seq_data[4:12, "Count"] <- 10
    seq_data[4:12, "FractionOfTotal"] <- 10 / tot
    seq_data[4:12, "FractionOfLocus"] <- 10 / tot
    sample_data <- analyze_sample(seq_data,  list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data,  list(Locus = "A"),
                                       counts.min = 500)
    with(sample_summary, {
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
      expect_equal(Ambiguous, FALSE)
      expect_equal(Stutter, TRUE)
      expect_equal(Artifact, FALSE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, 4910)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample handles ambiguous sequences", {
    # If sequences contain unexpected characters these should be interpreted as
    # ambiguous and filtered.
    # Replace the second-highest-count sequence to include an ambiguous base
    # near the end.
    idx <- nchar(seqs3$A) == 170
    seqs3$A[idx] <- sub("AGCCAGTC$",
                        "AGCCNAGTC",
                        seqs3$A[idx])
    seq_data <- analyze_seqs(seqs3$A, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data, counts.min = 500)
    with(sample_summary, {
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
      expect_equal(Ambiguous, TRUE)
      expect_equal(Stutter, FALSE)
      expect_equal(Artifact, FALSE)
      expect_equal(CountTotal, 5000)
      expect_equal(CountLocus, 4500)
      expect_equal(ProminentSeqs, 1)
    })
  })

  test_that("summarize_sample counts prominent sequences", {
    seq_data_1B          <- analyze_seqs(seqs1$B, locus_attrs, 3)
    sample_data_1B       <- analyze_sample(seq_data_1B, list(Locus = "B"), 0.05)
    sample_summary_1B    <- summarize_sample(sample_data_1B,
                                             list(Locus = "B"),
                                             counts.min = 500)
    seq_data_2B          <- analyze_seqs(seqs2$B, locus_attrs, 3)
    sample_data_2B       <- analyze_sample(seq_data_2B, list(Locus = "B"), 0.05)
    sample_summary_2B    <- summarize_sample(sample_data_2B,
                                             list(Locus = "B"),
                                             counts.min = 500)
    seq_data_3B          <- analyze_seqs(seqs3$B, locus_attrs, 3)
    sample_data_3B       <- analyze_sample(seq_data_3B, list(Locus = "B"), 0.05)
    sample_summary_3B    <- summarize_sample(sample_data_3B,
                                             list(Locus = "B"),
                                             counts.min = 500)
    seq_data_empty       <- analyze_seqs(c(), locus_attrs, 3)
    sample_data_empty    <- analyze_sample(seq_data_empty,
                                           list(Locus = "B"),
                                           0.05)
    sample_summary_empty <- summarize_sample(sample_data_empty,
                                             list(Locus = "B"),
                                             counts.min = 500)

    expect_equal(sample_summary_1B$ProminentSeqs,    2)
    expect_equal(sample_summary_2B$ProminentSeqs,    3)
    expect_equal(sample_summary_3B$ProminentSeqs,    1)
    expect_equal(sample_summary_empty$ProminentSeqs, 0)
    # Despite having stutter-y peaks the first two did not have a potential
    # allele removed, so Stutter == FALSE.  The third had stutter large enough
    # to look like an allele so it was removed.  The last one was empty so there
    # was nothing to have stutter from.
    expect_equal(sample_summary_1B$Stutter,    FALSE)
    expect_equal(sample_summary_2B$Stutter,    FALSE)
    expect_equal(sample_summary_3B$Stutter,    TRUE)
    expect_equal(sample_summary_empty$Stutter, FALSE)
  })

  test_that("summarize_sample rejects low-count samples", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    # Here we check that the filtered-counts-thresholding is applied, by forcing
    # the counts to a low number.  This should still report some stats but
    # should leave out the allele1/allele2 information.
    seq_data$Count <- seq_data$Count / 100
    sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data, list(Locus = "A"),
                                       counts.min = 500)
    with(sample_summary, {
      expect_equal(Allele1Seq, as.character(NA))
      expect_equal(Allele1Count, as.integer(NA))
      expect_equal(Allele1Length, as.integer(NA))
      expect_equal(Allele2Seq, as.character(NA))
      expect_equal(Allele2Count, as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous, FALSE)
      expect_equal(Ambiguous, FALSE)
      expect_equal(Stutter, FALSE)
      expect_equal(Artifact, FALSE)
      expect_equal(CountTotal, 50)
      expect_equal(CountLocus, 45)
      expect_equal(ProminentSeqs, 2)
    })
  })

  test_that("summarize_sample works with vector for sample attrs", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, c(Locus = "A"), 0.05)
    sample_summary <- summarize_sample(sample_data, c(Locus = "A"),
                                       counts.min = 500)
    check.seqs1A_summary(sample_summary)
  })

  test_that("summarize_sample warns of missing locus name", {
    # summarize_sample() should be able to tell if an invalid (as per the
    # earlier processing) locus name is given, because it won't be in the levels
    # of the MatchingLocus factor of the data frame from analyze_seqs().  I
    # think this would only come about when locus names given by
    # prepare.dataset() don't match what's in locus_attrs.
    skip("test not yet implemented")
  })

# test summarize_sample_guided --------------------------------------------

  # So do the summarize_ functions take lists or vectors for sample.attrs?  Make
  # that explicit and consistent!

  test_that("summarize_sample_guided summarizes sample attributes", {
    # This should have the same behavior as summarize_sample above.
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample_guided(sample_data, list(Locus = "A"),
                                      counts.min = 500)
    check.seqs1A_summary(sample_summary)
  })

  test_that("summarize_sample_guided handles completely empty sample data", {
    # This should have the same behavior as summarize_sample above.
    seq_data <- analyze_seqs(c(), locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, list(Locus = "A"), 0.05)
    sample_summary <- summarize_sample_guided(sample_data,  list(Locus = "A"),
                                       counts.min = 500)
    expect_equal(names(sample_summary), sample.summary.cols)
    expect_empty(sample_summary)
  })

  test_that("summarize_sample_guided uses expected_lengths", {
    # This should give allele sequences matching the given expected_lengths,
    # including order.
    # Flip the order of alleles here to check that aspect
    sample.attrs <- list(Locus = "A",
                         ExpectedLength1 = 194,
                         ExpectedLength2 = 162)
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 500)
    check.seqs1A_summary(sample_summary, ord = 2:1)
  })

  test_that("summarize_sample_guided ignores counts.min for expected lengths", {
    # This should give allele sequences matching the given expected_lengths,
    # including order, despite the counts.min value.
    # Flip the order of alleles here to check that aspect
    sample.attrs <- list(Locus = "A",
                         ExpectedLength1 = 194,
                         ExpectedLength2 = 162)
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 5000)
    check.seqs1A_summary(sample_summary, ord = 2:1)
  })

  test_that("summarize_sample_guided ignores counts.min for one exp. length", {
    # This should give allele sequences matching the given expected_lengths,
    # including order, despite the counts.min value.
    sample.attrs <- list(Locus = "A",
                         ExpectedLength1 = 194,
                         ExpectedLength2 = 194)
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 5000)
    # Check for one and only one called allele.
    with(sample_summary, {
      expect_equal(Allele1Seq,
                   gsub("[\n ]*", "",
                        "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAG
                        ATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
                        TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT
                        AGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC"))
      expect_equal(Allele2Seq,    as.character(NA))
      expect_equal(Allele1Count,  1300)
      expect_equal(Allele1Length, 194)
      expect_equal(Allele2Count,  as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous,    TRUE)
      expect_equal(Ambiguous,     FALSE)
      expect_equal(Stutter,       FALSE)
      expect_equal(Artifact,      FALSE)
      expect_equal(CountTotal,    5000)
      expect_equal(CountLocus,    4500)
      expect_equal(ProminentSeqs, 1)
    })

  })

  test_that("summarize_sample_guided uses counts.min if no expected lengths", {
    # This should not report alleles if total filtered read count is below
    # counts.min threshold.
    sample.attrs <- list(Locus = "A",
                         ExpectedLength1 = NA,
                         ExpectedLength2 = NA)
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 5000)
    with(sample_summary, {
      expect_equal(Allele1Seq,    as.character(NA))
      expect_equal(Allele1Count,  as.integer(NA))
      expect_equal(Allele1Length, as.integer(NA))
      expect_equal(Allele2Seq,    as.character(NA))
      expect_equal(Allele2Count,  as.integer(NA))
      expect_equal(Allele2Length, as.integer(NA))
      expect_equal(Homozygous,    FALSE)
      expect_equal(Ambiguous,     FALSE)
      expect_equal(Stutter,       FALSE)
      expect_equal(Artifact,      FALSE)
      expect_equal(CountTotal,    5000)
      expect_equal(CountLocus,    4500)
      expect_equal(ProminentSeqs, 2)
    })
  })

  test_that("summarize_sample_guided works with vector for sample attrs", {
    # It's "supposed" to be a list, but when used with functions like apply with
    # a data frame it gets munged into a vector.  So let's make sure that works
    # too.
    # Flip the order of alleles here to check that aspect too
    sample.attrs <- unlist(list(Locus = "A",
                                ExpectedLength1 = 194,
                                ExpectedLength2 = 162))
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 500)
    check.seqs1A_summary(sample_summary, ord = 2:1)
  })

  test_that("summarize_sample_guided works with blank expected_lengths", {
    # This should have the same behavior as summarize_sample above.
    sample.attrs <- list(Locus = "A",
                         ExpectedLength1 = NA,
                         ExpectedLength2 = NA)
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    sample_data <- analyze_sample_guided(seq_data, sample.attrs, 0.05)
    sample_summary <- summarize_sample_guided(sample_data, sample.attrs,
                                              counts.min = 500)
    check.seqs1A_summary(sample_summary)
  })

})
