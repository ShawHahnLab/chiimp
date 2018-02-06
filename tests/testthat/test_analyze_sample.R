context("Test sample analysis")

with(test_data, {

# test analyze_sample -----------------------------------------------------

  test_that("analyze_sample tabulates sequences", {
    sample.data <- analyze_sample(seqs1$A, locus_attrs, 3)
    expect_equal(nrow(sample.data), 14)
    chunk <- sample.data[1:2, ]
    expect_equal(chunk[1, "Length"], 162)
    expect_equal(chunk[2, "Length"], 194)
    with(chunk, {
      expect_equal(droplevels(MatchingLocus), factor(c("A", "A")))
      expect_equal(MotifMatch, c(T, T))
      expect_equal(LengthMatch, c(T, T))
      expect_equal(Stutter, as.integer(c(NA, NA)))
    })
  })

  test_that("analyze_sample handles completely emtpy input vector", {
    sample.data <- analyze_sample(c(), locus_attrs, 3)
    expect_equal(nrow(sample.data), 0)
    expect_equal(colnames(sample.data), sample.data.cols)
  })

  test_that("analyze_sample handles empty sequences", {
    seqs <- seqs1$A
    seqs[1:100] <- "" # empty out a segment of the vector
    sample.data <- analyze_sample(seqs, locus_attrs, 3)
    # Looking specifically at the entry for zero length
    chunk <- subset(sample.data, Length == 0)
    # There should be one row, accounting for all 100 original blank entries
    expect_equal(nrow(chunk), 1)
    expect_equal(chunk$Seq, "")
    expect_equal(chunk$Count, 100)
    # There should be nothing in the other columns
    expect_equal(droplevels(chunk$MatchingLocus), factor(NA, levels = c()))
    expect_equal(chunk$MotifMatch, FALSE)
    expect_equal(chunk$LengthMatch, NA)
    expect_equal(chunk$Stutter, as.integer(NA))
  })

  test_that("analyze_sample checks for motif repeats", {
    seqs <- seqs1$A
    sample.data <- analyze_sample(seqs, locus_attrs, 3)
    chunk <- subset(sample.data, !MotifMatch)
    with(chunk, {
      expect_equal(sum(Count), 500)
      expect_equal(range(Length), c(45, 57))
      expect_equal(droplevels(unique(MatchingLocus)), factor("A"))
      expect_true(all(is.na(Stutter)))
    })
  })

  test_that("analyze_sample checks for length", {
    seqs <- seqs1$A
    sample.data <- analyze_sample(seqs, locus_attrs, 3)
    chunk <- subset(sample.data, !LengthMatch)
    with(chunk, {
    expect_equal(sum(Count), 500)
    expect_equal(range(Length), c(45, 57))
    expect_equal(droplevels(unique(MatchingLocus)), factor("A"))
    expect_true(all(!(MotifMatch)))
    expect_true(all(is.na(Stutter)))
    })
  })

  test_that("analyze_sample marks stutter entries", {
    sample.data <- analyze_sample(seqs1$A, locus_attrs, 3)
    chunk <- subset(sample.data, !is.na(Stutter))
    expect_equal(chunk$Count, c(281, 116))
    expect_equal(chunk$Stutter, c(1, 2))
    # TODO add expectation for high-count case that should not be marked as
    # stutter
  })

})
