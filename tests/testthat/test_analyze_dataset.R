context("Test dataset analysis")

# Cases to cover: fastq, fasta, fasta.gz, fastq.gz, multiplex, ...?
# At this point we should have all the basics covered, and can cover the further
# analysis in summarize_dataset.

# Currently if I try to run analyze_dataset with the parallel package the tests
# fail when run under R CMD check but work fine when run interactively.  For now
# I'm just sticking with ncores = 1 here and in test_summarize_dataset.R to
# avoid calling parallel:: functions.
# Possibly relevant:
#  * https://github.com/r-lib/testthat/issues/602
#  * https://github.com/hadley/devtools/issues/1526
#  * https://github.com/r-lib/testthat/issues/86

with(test_data, {

# test analyze_dataset ----------------------------------------------------

  test_that("analyze_dataset produces expected list structure", {
    # Preliminary check on the data structure returned.
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    # prepare_dataset tested separately in test_io.R
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    results <- analyze_dataset(dataset, locus_attrs,
                               analysis_opts = list(fraction.min = 0.05),
                               summary_opts = list(counts.min = 500),
                               nrepeats = 3,
                               ncores = 1)
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)
    # Check the overall structure
    expect_equal(sapply(results, class),
                 c(summary = "data.frame",
                   samples = "list",
                   files = "list"))
  })

  test_that("analyze_dataset processes samples correctly", {
    # The general case for analyze_dataset.
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    results <- analyze_dataset(dataset, locus_attrs,
                               analysis_opts = list(fraction.min = 0.05),
                               summary_opts = list(counts.min = 500),
                               nrepeats = 3,
                               ncores = 1)
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)
    # Check the summary data frame
    with(results$summary, {
      # First update ordering of dataset's rows.  The existing order should be
      # correct except for the locus order definied via locus_attrs.
      dataset <- dataset[order(match(dataset$Locus, locus_attrs$Locus)), ]
      # Now, these should just match what we fed in via dataset.
      expect_equal(Filename,  dataset$Filename)
      expect_equal(Replicate, dataset$Replicate)
      expect_equal(Sample,    dataset$Sample)
      expect_equal(as.character(Locus), dataset$Locus)
      # Skipping allele seqs here, but we should already be checking more than
      # enough and have that checked in test_summarize_sample.
      expect_equal(Allele1Count, c(2783, 2261, 3415, 2183, 2703, 2336,
                                   2809, 2770, 2754, 2117, 4136, 2615))
      expect_equal(Allele1Length, c(162, 162, 182, 252, 236, 236,
                                    280, 284, 280, 250, 266, 342))
      expect_equal(Allele2Count, c(1290, 1309, NA, 967, 890, 974,
                                   1078, 1047, 1382, 1146, NA, 1217))
      expect_equal(Allele2Length, c(194, 178, NA, 216, 240, 220,
                                    260, 256, 276, 318, NA, 238))
      # Auto-generated names for all sequences.  NA entries in give NA out.
      expect_equal(Allele1Name, c("162-c6933c", "162-c6933c", "182-d679e1",
                                  "252-27c5bf", "236-321c79", "236-321c79",
                                  "280-74dd46", "284-2b3fab", "280-74dd46",
                                  "250-5dacee", "266-2aa675", "342-2e88c0"))
      expect_equal(Allele2Name, c("194-fc013a", "178-d84dc0", NA,
                                  "216-c0f11a", "240-2a344f", "220-fb9a92",
                                  "260-9a01fc", "256-c18a06", "276-ea279a",
                                  "318-35b7b6", NA, "238-6cc8ff"))
      h <- logical(12)
      h[c(3, 11)] <- TRUE
      expect_equal(Homozygous, h)
      s <- logical(12)
      s[c(3, 11)] <- TRUE
      expect_equal(Stutter, s)
      a <- logical(12)
      expect_equal(Artifact, a)
      expect_equal(CountTotal, integer(12) + 5000)
      expect_equal(CountLocus, integer(12) + 4466)
      expect_equal(ProminentSeqs,  c(2, 3, 1, 4, 2, 3, 2, 3, 2, 3, 1, 2))
    })

  })

  test_that("analyze_dataset names known alleles", {
    # If we gave names for some known allele sequences, do they show up
    # appropriately in the summary table?

    # First, set up example as above, but using known_alleles data frame
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    known_alleles <- data.frame(Locus = c("1", "1", "A"),
      Seq = c(paste0("ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTA",
                     "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
                     "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
                     "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
                     "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTGTGGCTCA",
                     "AAAGCTGAAT"),
              paste0("ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTA",
                     "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
                     "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
                     "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
                     "TATCTATCTATCTATCTATCTATCCTGTGGCTCAAAAGCTGAAT"),
              paste0("TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAGATAGA",
                     "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA",
                     "GATAGATAGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC")),
      Name = c("280-a", "260-X", "different_name_format"))
    results <- analyze_dataset(dataset, locus_attrs,
                               analysis_opts = list(fraction.min = 0.05),
                               summary_opts = list(counts.min = 500),
                               nrepeats = 3,
                               ncores = 1,
                               known_alleles = known_alleles)
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)

    # Check that the resulting allele names match all the expected values
    with(results$summary, {
      expect_equal(Allele1Name, c("different_name_format",
                                  "different_name_format",
                                  "182-d679e1", "252-27c5bf", "236-321c79",
                                  "236-321c79", "280-a", "284-2b3fab",
                                  "280-a", "250-5dacee", "266-2aa675",
                                  "342-2e88c0"))
      expect_equal(Allele2Name, c("194-fc013a", "178-d84dc0", NA,
                                  "216-c0f11a", "240-2a344f", "220-fb9a92",
                                  "260-X", "256-c18a06", "276-ea279a",
                                  "318-35b7b6", NA, "238-6cc8ff"))
    })

    # Check that the resulting sequence names in the sample data frames match
    # up. We have the one or two called alleles per sample checked here, with
    # one cross-appearance checked below.
    lapply(rownames(results$summary), function(nm) {
      # First called allele for these cases should always be the first seq in
      # each table.
      expect_equal(results$summary[nm, "Allele1Name"],
                   results$samples[[nm]]$SeqName[1])
      # Second called allele, if present, will be below the first somewhere.
      # Remaining seqs will be unnamed.
      if (! results$summary[nm, "Homozygous"]) {
        idx <- match(results$summary[nm, "Allele2Seq"],
                     results$samples[[nm]]$Seq)
        expect_equal(results$summary[nm, "Allele2Name"],
                     results$samples[[nm]]$SeqName[idx])
      }
    })

    # One particular case: 3-B showed a stutter-rejected sequence that's the
    # called allele for another sample.
    expect_equal(results$samples[["3-B"]]$SeqName[2], "220-fb9a92")
  })

  test_that("analyze_dataset handles missing loci", {
    # If there are locus names in dataset$Locus that are not present in the
    # rownames of locus_attrs, it should throw an error.
    data.dir <- tempfile()
    # the names are case-sensitive!
    seqs <- lapply(seqs, function(s) {
      names(s) <- c("a", "b", 1, 2)
      s
      })
    write_seqs(seqs, data.dir)
    # prepare_dataset tested separately in test_io.R
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    expect_error({
      results <- analyze_dataset(dataset, locus_attrs,
                                 analysis_opts = list(fraction.min = 0.05),
                                 summary_opts = list(counts.min = 500),
                                 nrepeats = 3,
                                 ncores = 1)
    }, "ERROR: Locus names in dataset not in attributes table: a, b")
  })

  test_that("analyze_dataset warns of empty input files", {
    # If we have no reads at all right from the start, we should warn the user.
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    # empty out one file
    fps <- list.files(data.dir, full.names = TRUE)
    unlink(fps[1])
    touch(fps[1])
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    msg <- capture.output({
      results <- analyze_dataset(dataset, locus_attrs,
                                 analysis_opts = list(fraction.min = 0.05),
                                 summary_opts = list(counts.min = 500),
                                 nrepeats = 3,
                                 ncores = 1)
    }, type = "message")
    msg_exp <- "WARNING: Zero reads for 1 of 12 data files"
    expect_true(length(grep(msg_exp, msg)) == 1)
  })

})
