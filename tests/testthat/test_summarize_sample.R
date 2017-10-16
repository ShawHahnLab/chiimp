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
  # would otherwise get counted (the rest were already set to off-target
  # junk during setup).
  check.seqs1A_summary(sample.summary,
                       count.locus = 4414,
                       allele1.count = 2748,
                       allele2.count = 1276)
})

test_that("summarize.sample marks stutter removal", {
  sample.data <- analyze.sample(seqs3$A, locus_attrs, 3)
  sample.summary <- summarize.sample(sample.data, "A")
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

test_that("summarize.sample counts prominent sequences", {
  sample.data.1B       <- analyze.sample(seqs1$B, locus_attrs, 3)
  sample.summary.1B    <- summarize.sample(sample.data.1B, "B")
  sample.data.2B       <- analyze.sample(seqs2$B, locus_attrs, 3)
  sample.summary.2B    <- summarize.sample(sample.data.2B, "B")
  sample.data.3B       <- analyze.sample(seqs3$B, locus_attrs, 3)
  sample.summary.3B    <- summarize.sample(sample.data.3B, 'B')
  sample.data.empty    <- analyze.sample(c(), locus_attrs, 3)
  sample.summary.empty <- summarize.sample(sample.data.empty, "B")

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

test_that("summarize.sample rejects low-count samples", {
  sample.data <- analyze.sample(seqs1$A, locus_attrs, 3)
  # Here we check that the filtered-counts-thresholding is applied, by forcing
  # the counts to a low number.  This should still report some stats but should
  # leave out the allele1/allele2 information.
  sample.data$Count <- sample.data$Count/100
  sample.summary <- summarize.sample(sample.data, "A")
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
