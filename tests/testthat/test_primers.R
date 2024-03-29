# make_read_primer_table --------------------------------------------------


test_that("make_read_primer_table finds locus matches", {
  reads <- c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTCATGGCATTGGCCGC",
    "AAAATGAACCCTCTAGCGACTGACCGTAGAACTCTGATCAGAAGTCAATCT",
    "TTCCACGCTCTGGGGAAAGGAATACGAAGTCACCTTGCTCCGTGTCACGGCTGATG",
    "ATCCAAGGGCATTACCACCCTAGCTCAATCGTCGCCGATTACATAG",
    "GGGGTTCTGCATGACGGGCCCGGAGCTATATGCGTTGGGGAGTGAACGC")
  # match three of the five reads
  primers <- c(
    "GGTCAGGTCAATTATTA",
    "CCACGCTCTGGGGAAAGGAA",
    "GGGGTTCTGCCTGACGG")
  # likewise, though not all the same ones
  primers_rev <- c(
    "CATGGCACTGGCCGC",
    "TCTGAAGTCAATCT",
    "ATGCGTTGGGGTG")
  locus_attrs <- data.frame(
    Locus = c("X", "Y", "Z"),
    Primer = primers,
    ReversePrimer = primers_rev,
    stringsAsFactors = FALSE)
  # should call match_primer_set for both fwd and rev, and then only keep those
  # with the same locus found for both fwd and reverse.
  result <- make_read_primer_table(
    reads, locus_attrs, max_mismatches = 3, primer_action = "none",
    use_reverse_primers = TRUE, reverse_primer_r1 = TRUE)
  expect_identical(
    colnames(result),
    c("SeqOrig",
      "FwdStart", "FwdStop", "FwdMismatches", "FwdLocus",
      "RevStart", "RevStop", "RevMismatches", "RevLocus",
      "MatchingLocus", "Seq"))
  expect_identical(result$SeqOrig, reads)
  expect_identical(result$MatchingLocus, c("X", NA, NA, NA, "Z"))
  expect_identical(result$Seq, reads)
})

test_that("make_read_primer_table can handle missing primer seqs", {
  reads <- c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTCATGGCATTGGCCGC",
    "AAAATGAACCCTCTAGCGACTGACCGTAGAACTCTGATCAGAAGTCAATCT",
    "TTCCACGCTCTGGGGAAAGGAATACGAAGTCACCTTGCTCCGTGTCACGGCTGATG",
    "ATCCAAGGGCATTACCACCCTAGCTCAATCGTCGCCGATTACATAG",
    "GGGGTTCTGCATGACGGGCCCGGAGCTATATGCGTTGGGGAGTGAACGC")
  primers <- c(
    "GGTCAGGTCAATTATTA",
    "CCACGCTCTGGGGAAAGGAA",
    "GGGGTTCTGCCTGACGG")
  primers_rev <- c(
    "CATGGCACTGGCCGC",
    "TCTGAAGTCAATCT",
    "")
  locus_attrs <- data.frame(
    Locus = c("X", "Y", "Z"),
    Primer = primers,
    ReversePrimer = primers_rev,
    stringsAsFactors = FALSE)
  result <- make_read_primer_table(
    reads, locus_attrs, max_mismatches = 3, primer_action = "none",
    use_reverse_primers = TRUE, reverse_primer_r1 = TRUE)
  expect_identical(result$SeqOrig, reads)
  expect_identical(result$MatchingLocus, c("X", NA, NA, NA, NA))
  expect_identical(result$Seq, reads)
  # Should treat NA like ""
  locus_attrs$ReversePrimer[3] <- NA
  result2 <-  make_read_primer_table(
    reads, locus_attrs, max_mismatches = 3, primer_action = "none",
    use_reverse_primers = TRUE, reverse_primer_r1 = TRUE)
  expect_identical(result2, result)
})

test_that("make_read_primer_table understands IUPAC codes in primer seqs", {
  reads <- c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CAGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CTGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CCGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC")
  locus_attrs <- data.frame(
    Locus = "X",
    Primer = "CNGCCAGTGGAGCGGGC", # perfect match at start of all reads
    stringsAsFactors = FALSE)
  # basic case, just an N in the second position that should match any base
  result <- make_read_primer_table(
    reads, locus_attrs, max_mismatches = 3,
    use_reverse_primers = FALSE)
  expect_identical(result$FwdStart, rep(1L, 4))
  expect_identical(result$FwdStop, rep(17L, 4))
  expect_identical(result$FwdMismatches, rep(0L, 4))
  expect_identical(result$FwdLocus, rep("X", 4))
  # but ambiguous bases are only allowed to match for primers, not reads
  result <- make_read_primer_table(
    "CNGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC", locus_attrs, 3,
    use_reverse_primers = FALSE)
  expect_identical(result$FwdStart, 1L)
  expect_identical(result$FwdStop, 17L)
  expect_identical(result$FwdMismatches, 1L)
  expect_identical(result$FwdLocus, "X")
})


# handle_primers ----------------------------------------------------------


test_that("handle_primers modifies sequences as directed", {
  # four cases for both fwd and rev:
  # none, keep, replace, remove
  # The primers don't exactly match what's in the sequences (maybe mismatches,
  # or IUPAC) but that shouldn't matter
  locus_attrs <- data.frame(
    Locus = "A",
    Primer = "XQXXX",
    ReversePrimer = "ZZZWZ",
    stringsAsFactors = FALSE)
  input <- data.frame(
    SeqOrig = "XXXXXYYYYYZZZZZ",
    FwdStart = 1,
    FwdStop = 5,
    RevStart = 11,
    RevStop = 15,
    MatchingLocus = "A",
    stringsAsFactors = FALSE)
  # general checks
  result <- handle_primers(input, locus_attrs, "none", "none", TRUE)
  expect_identical(
    result, cbind(input, Seq = input$SeqOrig, stringsAsFactors = FALSE))
  result <- handle_primers(input, locus_attrs, "keep", "keep", TRUE)
  expect_identical(
    result, cbind(input, Seq = input$SeqOrig, stringsAsFactors = FALSE))
  result <- handle_primers(input, locus_attrs, "replace", "remove", TRUE)
  expect_identical(
    result, cbind(input, Seq = "XQXXXYYYYY", stringsAsFactors = FALSE))
  result <- handle_primers(input, locus_attrs, "remove", "replace", TRUE)
  expect_identical(
    result, cbind(input, Seq = "YYYYYZZZWZ", stringsAsFactors = FALSE))
  result <- handle_primers(input, locus_attrs, "remove", "replace", FALSE)
  expect_identical(
    result, cbind(input, Seq = "YYYYYZWZZZ", stringsAsFactors = FALSE))
  # check trimming; if primers are set in a bit, that should show in the output
  input$FwdStart <- 2
  input$RevStop <- 12
  result <- handle_primers(input, locus_attrs, "none", "none", TRUE)
  expect_identical(
    result, cbind(input, Seq = input$SeqOrig, stringsAsFactors = FALSE))
  result <- handle_primers(input, locus_attrs, "keep", "keep", TRUE)
  expect_identical(
    result, cbind(input, Seq = "XXXXYYYYYZZ", stringsAsFactors = FALSE))
})

test_that("handle_primers handles empty inputs", {
  locus_attrs <- data.frame(
    Locus = "A",
    Primer = "XQXXX",
    ReversePrimer = "ZZZWZ",
    stringsAsFactors = FALSE)
  input <- data.frame(
    SeqOrig = "",
    FwdStart = NA,
    FwdStop = NA,
    RevStart = NA,
    RevStop = NA,
    MatchingLocus = NA,
    stringsAsFactors = FALSE)
  # empty strings should be fine
  result <- handle_primers(input, locus_attrs, "none", "none", TRUE)
  expected <- cbind(input, Seq = input$SeqOrig, stringsAsFactors = FALSE)
  expect_identical(result, expected)
  # or NA
  input$SeqOrig <- as.character(NA)
  expected <- cbind(input, Seq = input$SeqOrig, stringsAsFactors = FALSE)
  result <- handle_primers(input, locus_attrs, "none", "none", TRUE)
  expect_identical(result, expected)
  # no rows in, no rows out
  result <- handle_primers(input[0, ], locus_attrs, "none", "none", TRUE)
  expected <- cbind(
    input[0, ], Seq = input$SeqOrig[0], stringsAsFactors = FALSE)
  expect_identical(result, expected)
})

test_that("handle_primers handles invalid keywords", {
  locus_attrs <- data.frame(
    Locus = "A",
    Primer = "XQXXX",
    ReversePrimer = "ZZZWZ",
    stringsAsFactors = FALSE)
  input <- data.frame(
    SeqOrig = "XXXXXYYYYYZZZZZ",
    FwdStart = 1,
    FwdStop = 5,
    RevStart = 11,
    RevStop = 15,
    MatchingLocus = "A",
    stringsAsFactors = FALSE)
  expect_error(
    handle_primers(input, locus_attrs, "trim", "none", TRUE),
    "primer_action_fwd should be")
  expect_error(
    handle_primers(input, locus_attrs, "none", "trim", TRUE),
    "primer_action_rev should be")
})

test_that("handle_primers handles missing reverse primer info", {
  # If make_read_primer_table is called with use_reverse_primers=FALSE, all the
  # Rev* columns will be missing.  handle_primers will need to handle this case.
  locus_attrs <- data.frame(
    Locus = "A",
    Primer = "XQXXX",
    ReversePrimer = "ZZZWZ",
    stringsAsFactors = FALSE)
  input <- data.frame(
    SeqOrig = "XXXXXYYYYY",
    FwdStart = 1,
    FwdStop = 5,
    MatchingLocus = "A",
    stringsAsFactors = FALSE)
  # this should halt and complain
  expect_error(
    handle_primers(input, locus_attrs, "remove", "remove"),
    "can't apply primer_action_rev")
  # but this should work
  result <- handle_primers(input, locus_attrs, "remove", "none")
  expected <- cbind(input, Seq = "YYYYY", stringsAsFactors = FALSE)
  expect_identical(result, expected)
})

# find_primer_matches -----------------------------------------------------


test_that("find_primer_matches tabulates all the best primer matches", {
  # The result should be all reads compared to all primers, with SeqIdx
  # corresponding to each read.
  reads <- c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTCATGGCATTGGCCGC",
    "AAAATGAACCCTCTAGCGACTGACCGTAGAACTCTGATCAGAAGTCAATCT",
    "TTCCACGCTCTGGGGAAAGGAATACGAAGTCACCTTGCTCCGTGTCACGGCTGATG",
    "ATCCAAGGGCATTACCACCCTAGCTCAATCGTCGCCGATTACATAG",
    "GGGGTTCTGCATGACGGGCCCGGAGCTATATGCGTTGGGGAGTGAACGC")
  primers <- c(
    "CACGTTAGAGCTCATTGT",
    "ACCGTCAAAGCGCTCCGGGAGTACAGA",
    "ATAATACCGCGCCTAAGATTC")
  result <- find_primer_matches(reads, primers, NA)
  expect_identical(
    colnames(result),
    c("SeqIdx", "PrimerIdx", "Start", "Stop", "Mismatches"))
  expect_identical(result$SeqIdx, rep(seq_along(reads), times = 3))
  expect_identical(result$PrimerIdx, rep(seq_along(primers), each = 5))
   expect_identical(
     result$Mismatches,
     as.integer(c(9, 8, 9, 10, 11, 17, 16, 16, 15, 15, 13, 12, 13, 12, 12)))
})

test_that("find_primer_matches handles empty inputs", {
  # If one or the other (or both) are zero-length, the result should be
  # an empty data frame
  reads <- c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTCATGGCATTGGCCGC",
    "AAAATGAACCCTCTAGCGACTGACCGTAGAACTCTGATCAGAAGTCAATCT",
    "TTCCACGCTCTGGGGAAAGGAATACGAAGTCACCTTGCTCCGTGTCACGGCTGATG",
    "ATCCAAGGGCATTACCACCCTAGCTCAATCGTCGCCGATTACATAG",
    "GGGGTTCTGCATGACGGGCCCGGAGCTATATGCGTTGGGGAGTGAACGC")
  primers <- c(
    "CACGTTAGAGCTCATTGT",
    "ACCGTCAAAGCGCTCCGGGAGTACAGA",
    "ATAATACCGCGCCTAAGATTC")
  results <- list(
    find_primer_matches(reads, character(), NA),
    find_primer_matches(character(), primers, NA),
    find_primer_matches(character(), character(), NA),
    # NAs and empty strings are excluded
    find_primer_matches(reads, substr(primers, 0, 0), NA),
    find_primer_matches(substr(reads, 0, 0), primers, NA),
    find_primer_matches(substr(reads, 0, 0), substr(primers, 0, 0), NA),
    find_primer_matches(reads, rep(as.character(NA), 3), NA),
    find_primer_matches(rep(as.character(NA), 5), primers, NA),
    find_primer_matches(rep(as.character(NA), 5), rep(as.character(NA), 3), NA))
  for (result in results) {
    expect_identical(
      colnames(result),
      c("SeqIdx", "PrimerIdx", "Start", "Stop", "Mismatches"))
    expect_equal(nrow(result), 0)
  }
})

test_that("find_primer_matches finds matches at edges", {
  # one match right at the start
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
    "TAAGAAA", NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1L,
      PrimerIdx = 1L,
      Start = 1L,
      Stop = 7L,
      Mismatches = 0L))
  # one match right at the end
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
                           "AATCAAC", NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1L,
      PrimerIdx = 1L,
      Start = 24L,
      Stop = 30L,
      Mismatches = 0L))
  # weirder case: one match right at end for one seq, for the shorter of two
  # available primers.  (Testing this because I initially screwed up the
  # boundaries for the sequence comparisons for multiple primers and want to
  # make sure it doesn't regress.)
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
    c(
      "AATCAAC",
      "GCTGTATCTA"), NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = c(1L, 1L),
      PrimerIdx = c(1L, 2L),
      Start = c(24L, 6L),
      Stop = c(30L, 15L),
      Mismatches = c(0L, 5L)))
  # even weirder: primer extends past the edge of the match by 1
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAA",
    c(
      "AATCAAC",
      "GCTGTATCTA"), NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = c(1L, 1L),
      PrimerIdx = c(1L, 2L),
      Start = c(24L, 6L),
      Stop = c(29L, 15L),
      Mismatches = c(1L, 5L)))
})

test_that("find_primer_matches can use maximum mismatch count", {
  # one match right at the start, searching exhaustively
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
    c("TAAGAAA", "ACTGTGA"), NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = c(1L, 1L),
      PrimerIdx = c(1L, 2L),
      Start = c(1L, 6L),
      Stop = c(7L, 12L),
      Mismatches = c(0L, 4L)))
  # one match right at the start, but with max_mismatches, it won't search
  # exhaustively and we just get the one exact match back
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
    c("TAAGAAA", "ACTGTGA"), max_mismatches = 1)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1L,
      PrimerIdx = 1L,
      Start = 1L,
      Stop = 7L,
      Mismatches = 0L))
  # 2 mismatches now, so it won't find a match at all
  result <- find_primer_matches(
    "TAAGAAATGCTTATATGGCCATAAATCAAC",
    c("TGGGAAA", "ACTGTGA"), max_mismatches = 1)
  expect_identical(
    result,
    data.frame(
      SeqIdx = integer(),
      PrimerIdx = integer(),
      Start = integer(),
      Stop = integer(),
      Mismatches = integer()))
})

test_that("find_primer_matches handles IUPAC codes", {
  # perfect match with the "N" in the second position of the primer
  result <- find_primer_matches(c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CAGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CTGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CCGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC"),
    "CNGCCAGTGGAGCGGGC", NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1:4,
      PrimerIdx = 1L,
      Start = 1L,
      Stop = 17L,
      Mismatches = 0L))
  # but ambiguous bases are only allowed to match for primers, not reads, so
  # this gets marked as a mismatch
  result <- find_primer_matches(
    "CNGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CNGCCAGTGGAGCGGGC", NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1L,
      PrimerIdx = 1L,
      Start = 1L,
      Stop = 17L,
      Mismatches = 1L))
  # a different IUPAC code so only some match
  result <- find_primer_matches(c(
    "CGGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CAGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CTGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC",
    "CCGCCAGTGGAGCGGGCAGGTCAATTATTAGAGTCTAGACTC"),
    "CMGCCAGTGGAGCGGGC", NA)
  expect_identical(
    result,
    data.frame(
      SeqIdx = 1:4,
      PrimerIdx = 1L,
      Start = 1L,
      Stop = 17L,
      Mismatches = c(1L, 0L, 1L, 0L)))
})


# make_raw_nt -------------------------------------------------------------


test_that("make_raw_nt makes raw matrix from text", {
  # three basic categories are recognized characters, unrecognized characters,
  # and pad positions (for shorter seqs)
  seqs <- c("ACTX", "", "ACG")
  out <- make_raw_nt(seqs)
  out_exp <- matrix(as.raw(c(
    0x01, 0x02, 0x08, 0x00,
    0x80, 0x80, 0x80, 0x80,
    0x01, 0x02, 0x04, 0x80)), ncol = 3)
  expect_identical(out, out_exp)
  # The pad and other values can be customized
  out <- make_raw_nt(seqs, pad = 0x40, other = 0x80)
  out_exp[out_exp == 0x80] <- as.raw(0x40)
  out_exp[out_exp == 0x00] <- as.raw(0x80)
  expect_identical(out, out_exp)
  # raw input for those should work too
  out <- make_raw_nt(seqs, pad = as.raw(0x40), other = as.raw(0x80))
  expect_identical(out, out_exp)
  # custom mapping should be supported
  map <- as.raw(1:4)
  names(map) <- c("T", "G", "C", "A")
  out <- make_raw_nt(seqs, map)
  out_exp <- matrix(as.raw(c(
    0x04, 0x03, 0x01, 0x00,
    0x80, 0x80, 0x80, 0x80,
    0x04, 0x03, 0x02, 0x80)), ncol = 3)
  expect_identical(out, out_exp)
})

test_that("make_raw_nt handles empty input", {
  # empty input, empty output:
  # one sequence (cols) with zero positions (rows)
  expect_equal(dim(make_raw_nt("")), c(0, 1))
  # zero-length input, zero-length output
  expect_equal(dim(make_raw_nt(character())), c(0, 0))
})

test_that("make_raw_nt handles IUPAC codes", {
  seqs <- c("SCNG", "", "ACG")
  out <- make_raw_nt(seqs)
  out_exp <- matrix(as.raw(c(
    6, 2, 15, 4,
    128, 128, 128, 128,
    1, 2, 4, 128)), ncol = 3)
  expect_identical(out, out_exp)
})
