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
    dataset <- prepare_dataset(data.dir, '()(\\d+)-([A-Za-z0-9]+).fasta')
    results <- analyze_dataset(dataset, locus_attrs,
                               summary_args = list(
                                 fraction.min = 0.05,
                                 counts.min = 500),
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
    dataset <- prepare_dataset(data.dir, '()(\\d+)-([A-Za-z0-9]+).fasta')
    results <- analyze_dataset(dataset, locus_attrs,
                               summary_args = list(
                                 fraction.min = 0.05,
                                 counts.min = 500),
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
      expect_equal(Allele1Count, c(2803, 2394, 3971, 2489, 2656, 4082,
                                   2794, 2041, 2909, 3035, 2902, 2288))
      expect_equal(Allele1Length, c(162, 182, 174, 244, 212, 224,
                                    272, 276, 276, 278, 346, 318))
      expect_equal(Allele2Count, c(1300, 2003,   NA, 1304, 1059,   NA,
                                   1055, 1402, 1004, 1142, 1088, 1270))
      expect_equal(Allele2Length, c(194, 178,  NA, 220, 220,  NA,
                                    260, 288, 252, 342, 314, 350))
      # Auto-generated names for all sequences.  NA entries in give NA out.
      expect_equal(Allele1Name, c("162-c6933c", "182-d679e1", "174-8a43ea",
                                  "244-3c2ff2", "212-6d4afb", "224-053ec2",
                                  "272-292a2a", "276-ea279a", "276-ea279a",
                                  "278-ae70f3", "346-b05233", "318-35b7b6"))
      expect_equal(Allele2Name, c("194-fc013a", "178-d84dc0", NA,
                                  "220-fb9a92", "220-fb9a92", NA,
                                  "260-9a01fc", "288-201179", "252-a5eee8",
                                  "342-2e88c0", "314-ce2338", "350-4acdbb"))
      h <- logical(12)
      h[c(3, 6)] <- TRUE
      expect_equal(Homozygous, h)
      s <- logical(12)
      s[c(3, 6)] <- TRUE
      expect_equal(Stutter, s)
      a <- logical(12)
      expect_equal(Artifact, a)
      expect_equal(CountTotal, integer(12)+5000)
      expect_equal(CountLocus, integer(12)+4500)
      expect_equal(ProminentSeqs,  c(2, 2, 1, 2, 3, 1, 2, 3, 2, 2, 2, 3))
    })

  })

  test_that("analyze_dataset names known alleles", {
    # If we gave names for some known allele sequences, do they show up
    # appropriately in the summary table?

    # First, set up example as above, but using known_alleles data frame
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    dataset <- prepare_dataset(data.dir, '()(\\d+)-([A-Za-z0-9]+).fasta')
    known_alleles <- data.frame(Locus = c("1", "1", "A"),
                                Seq = c("ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTGTGGCTCAAAAGCTGAAT",
                                        "ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTGTGGCTCAAAAGCTGAAT",
                                        "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC"),
                                Name = c("272-a",
                                         "260-X",
                                         "different_name_format"))
    results <- analyze_dataset(dataset, locus_attrs,
                               summary_args = list(
                                 fraction.min = 0.05,
                                 counts.min = 500),
                               nrepeats = 3,
                               ncores = 1,
                               known_alleles = known_alleles)
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)

    # Check that the resulting allele names match all the expected values
    with(results$summary, {
      # expect_equal(Allele1Name, c("different_name_format", "182-d679e1", "174-8a43ea",
      #                             NA, NA, NA,
      #                             "272-a", "276-ea279a", "276-ea279a",
      #                             "278-ae70f3", "346-b05233", "318-35b7b6"))
      expect_equal(Allele1Name, c("different_name_format", "182-d679e1", "174-8a43ea",
                                  "244-3c2ff2", "212-6d4afb", "224-053ec2",
                                  "272-a", "276-ea279a", "276-ea279a",
                                  "278-ae70f3", "346-b05233", "318-35b7b6"))
      expect_equal(Allele2Name, c("194-fc013a", "178-d84dc0", NA,
                                  "220-fb9a92", "220-fb9a92", NA,
                                  "260-X", "288-201179", "252-a5eee8",
                                  "342-2e88c0", "314-ce2338", "350-4acdbb"))
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
        idx <- match(results$summary[nm, "Allele2Seq"], results$samples[[nm]]$Seq)
        expect_equal(results$summary[nm, "Allele2Name"],
                     results$samples[[nm]]$SeqName[idx])
      }
    })

    # One particular case: 3-B showed a stutter-rejected sequence that's the
    # called allele for another sample.
    expect_equal(results$samples[["3-B"]]$SeqName[2], "220-fb9a92")
  })

  test_that("analyze_dataset warns of missing loci", {
    # If there are locus names in dataset$Locus that are not present in the
    # rownames of locus_attrs, it should throw a warning.
    skip("test not yet implemented")
  })

})
