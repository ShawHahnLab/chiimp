context("Test sample analysis")

with(test_data, {


# test analyze_seqs -------------------------------------------------------


  test_that("analyze_seqs tabulates sequences", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    expect_equal(nrow(seq_data), 24)
    chunk <- seq_data[1:2, ]
    expect_equal(chunk[1, "Length"], 162)
    expect_equal(chunk[2, "Length"], 194)
    with(chunk, {
      expect_equal(droplevels(MatchingLocus), factor(c("A", "A")))
      expect_equal(MotifMatch, c(T, T))
      expect_equal(LengthMatch, c(T, T))
      expect_equal(Ambiguous, c(F, F))
      expect_equal(Stutter, as.integer(c(NA, NA)))
      expect_equal(Artifact, as.integer(c(NA, NA)))
    })
  })

  test_that("analyze_seqs handles completely emtpy input vector", {
    seq_data <- analyze_seqs(c(), locus_attrs, 3)
    expect_equal(nrow(seq_data), 0)
    expect_equal(colnames(seq_data), sample.data.cols)
  })

  test_that("analyze_seqs handles empty sequences", {
    seqs <- seqs1$A
    seqs[1:100] <- "" # empty out a segment of the vector
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    # Looking specifically at the entry for zero length
    chunk <- subset(seq_data, Length == 0)
    # There should be one row, accounting for all 100 original blank entries
    expect_equal(nrow(chunk), 1)
    expect_equal(chunk$Seq, "")
    expect_equal(chunk$Count, 100)
    # There should be nothing in the other columns
    expect_equal(droplevels(chunk$MatchingLocus), factor(NA, levels = c()))
    expect_equal(chunk$MotifMatch, FALSE)
    expect_equal(chunk$LengthMatch, NA)
    expect_equal(chunk$Ambiguous, FALSE)
    expect_equal(chunk$Stutter, as.integer(NA))
    expect_equal(chunk$Artifact, as.integer(NA))
  })

  test_that("analyze_seqs checks for motif repeats", {
    seqs <- seqs1$A
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    chunk <- subset(seq_data, !MotifMatch)
    with(chunk, {
      expect_equal(sum(Count), 483)
      expect_equal(range(Length), c(45, 57))
      expect_equal(droplevels(unique(MatchingLocus)), factor("A"))
      expect_true(all(is.na(Stutter)))
      expect_true(all(is.na(Artifact)))
    })
  })

  test_that("analyze_seqs checks for length", {
    seqs <- seqs1$A
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    chunk <- subset(seq_data, !LengthMatch)
    with(chunk, {
    expect_equal(sum(Count), 483)
    expect_equal(range(Length), c(45, 57))
    expect_equal(droplevels(unique(MatchingLocus)), factor("A"))
    expect_true(all(!(MotifMatch)))
    expect_true(all(is.na(Stutter)))
    expect_true(all(is.na(Artifact)))
    })
  })

  test_that("analyze_seqs marks stutter entries", {
    seq_data <- analyze_seqs(seqs1$A, locus_attrs, 3)
    chunk <- subset(seq_data, !is.na(Stutter))
    expect_equal(chunk$Count, c(279, 114, 2, 2, 2))
    expect_equal(chunk$Stutter, c(1, 2, 2, 2, 2))
    # TODO add expectation for high-count case that should not be marked as
    # stutter
  })

  test_that("analyze_seqs marks artifact entries", {
    s <- seqs1$A
    # Take that first stutter and make it an artifact instead
    highest <- names(sort(table(s), decreasing = T)[1])
    stutter <- names(sort(table(s), decreasing = T)[3])
    idx <- s == stutter
    s[idx] <- highest
    substr(s[idx], nchar(stutter), nchar(stutter)) <- "X"
    # Check that the third entry is marked an artifact of the first
    seq_data <- analyze_seqs(s, locus_attrs, 3)
    expect_equal(seq_data$Artifact, c(NA, NA, 1)[1:24])
  })

  test_that("analyze_seqs marks ambiguous entries", {
    # sequences that contain non-ACTG characters should be marked TRUE in the
    # Ambiguous column (interpreting those as "N" or similar).
    s <- seqs1$A
    s[s == s[1]] <- sub("AGCCAGTC", "AGCCANTC", s[1])
    seq_data <- analyze_seqs(s, locus_attrs, 3)
    expect_equal(seq_data$Ambiguous,
                 c(TRUE, rep(FALSE, nrow(seq_data) - 1)))
  })

})
