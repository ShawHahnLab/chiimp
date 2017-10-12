check.seqs1A_summary <- function(sample.summary, count.locus=4500) {
  expect_equal(names(sample.summary), sample.summary.cols)
  with(sample.summary, {
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
    expect_equal(Allele1Count, 2803)
    expect_equal(Allele1Length, 162)
    expect_equal(Allele2Count, 1300)
    expect_equal(Allele2Length, 194)
    expect_equal(Homozygous, FALSE)
    expect_equal(Stutter, FALSE)
    expect_equal(CountTotal, 5000)
    expect_equal(CountLocus, count.locus)
    expect_equal(ProminentSeqs, 3)
  })
}

test_that("summarize.sample summarizes sample attributes", {
  sample.data <- analyze.sample(seqs1$A, locus_attrs, 3)
  sample.summary <- summarize.sample(sample.data, "A")
  check.seqs1A_summary(sample.summary)
})

test_that("summarize.sample handles completely emtpy input sample data", {
  sample.data <- analyze.sample(c(), locus_attrs, 3)
  sample.summary <- summarize.sample(sample.data, "A")
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

test_that("summarize.sample handles empty sequences in input sample data", {
  seqs <- seqs1$A
  seqs[1:100] <- "" # empty out a segment of the vector
  sample.data <- analyze.sample(seqs, locus_attrs, 3)
  sample.summary <- summarize.sample(sample.data, "A")
  # Nothing should change in the output, except that we zeroed out 90 reads that
  # would otherwise get counted (the other 10 were already set to off-target
  # junk during setup).
  check.seqs1A_summary(sample.summary, count.locus = 4410)
})

test_that("summarize.sample marks stutter removal", {
  fail("test not yet implemented")
})

test_that("summarize.sample counts prominent sequences", {
  sample.data.1B       <- analyze.sample(seqs1$B, locus_attrs, 3)
  sample.summary.1B    <- summarize.sample(sample.data.1B, "B")
  sample.data.2B       <- analyze.sample(seqs2$B, locus_attrs, 3)
  sample.summary.2B    <- summarize.sample(sample.data.2B, "B")
  sample.data.3B       <- analyze.sample(seqs3$B, locus_attrs, 3)
  sample.summary.3B    <- summarize.sample(sample.data.3B, 'B')
  sample.data.empty    <- analyze.sample(c(), locus_attrs, 3)
  sample.summary.empty <- summarize.sample(sample.data.empty, "B")

  expect_equal(sample.summary.1B$ProminentSeqs,    1)
  expect_equal(sample.summary.2B$ProminentSeqs,    4)
  expect_equal(sample.summary.3B$ProminentSeqs,    3)
  expect_equal(sample.summary.empty$ProminentSeqs, 0)
  # TODO check cases where Stutter == TRUE
})

test_that("summarize.sample ignores low-count samples", {
  fail("feature not yet implemented")
  # Here we should check that the filtered-counts-thresholding (not yet
  # implemented!) is applied.
})
