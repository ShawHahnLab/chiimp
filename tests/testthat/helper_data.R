# simulated data for testing ----------------------------------------------


# I'm shoving all this into a list to keep it separate from all the other
# objects devtools::load_all() dumps into the namespace.  Users loading via
# devtools otherwise end up with a bunch of extra objects that are only relevant
# during testing.

#' Helper Data for Tests
#'
#' This list is a bundle of shared data and functions for running unit tests.
test_data <- within(list(), {
  txt.locus_attrs <- "Locus   LengthMin  LengthMax    LengthBuffer   Motif   Primer                  ReversePrimer
A        131        179          20             TAGA    TATCACTGGTGTTAGTCCTCTG  CACAGTTGTGTGAGCCAGTC
B        194        235          20             TAGA    AGTCTCTCTTTCTCCTTGCA    TAGGAGCCTGTGGTCCTGTT
1        232        270          20             TATC    ACAGTCAAGAATAACTGCCC    CTGTGGCTCAAAAGCTGAAT
2        218        337          20             TCCA    TTGTCTCCCCAGTTGCTA      TCTGTCATAAACCGTCTGCA"
  f.locus_attrs <- textConnection(txt.locus_attrs)
  locus_attrs <- read.table(f.locus_attrs, header = T, stringsAsFactors = F)
  rownames(locus_attrs) <- locus_attrs$Locus
  close(f.locus_attrs)
  rm(f.locus_attrs)

  sample.data.cols <- c("Seq", "Count", "Length", "MatchingLocus", "MotifMatch",
                        "LengthMatch", "Stutter", "Artifact", "FractionOfTotal",
                        "FractionOfLocus")
  sample.summary.cols <- c("Allele1Seq", "Allele1Count",
                           "Allele1Length", "Allele2Seq",
                           "Allele2Count", "Allele2Length",
                           "Homozygous", "Stutter", "Artifact", "CountTotal",
                           "CountLocus", "ProminentSeqs")

  make.seq_junk <- function(N) {
    nucleotides <- c("A", "T", "C", "G")
    vapply(runif(N, min = 1, max = 20), function(L)
      paste0(sample(nucleotides, L, replace = T), collapse = ""),
      "character")
  }

  simulate.seqs <- function(locus_name, locus_attrs, homozygous=NULL, N=5000,
                            off_target_ratio=10) {
    attrs <- locus_attrs[locus_name, ]
    L.min <- attrs$LengthMin - nchar(attrs$Primer)
    L.max <- attrs$LengthMax - nchar(attrs$Primer)
    L <- as.integer(runif(2, min = L.min, max = L.max))
    if (!missing(homozygous)) {
      if (homozygous) {
        L[2] <- L[1]
      } else {
        while (L[1] == L[2])
          L <- as.integer(runif(2, min = L.min, max = L.max))
      }
    }
    make.reps <- function(motif, len)
      paste(rep(motif, floor(len / nchar(motif))), collapse = "")
    repeats1 <- make.reps(attrs$Motif, L[1])
    repeats2 <- make.reps(attrs$Motif, L[2])
    seq.correct.1 <- paste0(attrs$Primer, repeats1, attrs$ReversePrimer)
    seq.correct.2 <- paste0(attrs$Primer, repeats2, attrs$ReversePrimer)
    # The exact correct sequences, either one or two unique
    N1 <- N / 2 + rnorm(10000) * N / 20
    N1 <- as.integer(max(min(N1, N * 0.9), N * 0.1))
    N2 <- N - N1
    seqs <- c(rep(seq.correct.1, N1),
              rep(seq.correct.2, N2))
    ## Mix in stutter
    N1.stutter <- as.integer(runif(1, min = N1 / 20, max = N1 / 3))
    N2.stutter <- as.integer(runif(1, min = N2 / 20, max = N2 / 3))
    seqs[1:N1.stutter] <- paste0(attrs$Primer,
                                 make.reps(attrs$Motif,
                                           L[1] - nchar(attrs$Motif)),
                                 attrs$ReversePrimer)
    seqs[N1 + 1:N2.stutter] <- paste0(attrs$Primer,
                                      make.reps(attrs$Motif,
                                                L[2] - nchar(attrs$Motif)),
                                      attrs$ReversePrimer)
    ## Mix in off target amplification
    idx <- seq(1, length(seqs), off_target_ratio)
    seqs[idx] <- paste0(attrs$Primer, make.seq_junk(10), attrs$ReversePrimer)
    ## TODO add contam
    ## TODO munge up the reads a bit

    ## Shuffle vector
    seqs <- sample(seqs)
    return(seqs)
  }

  simulate.set <- function(locus_attrs) {
    seqs <- lapply(rownames(locus_attrs), function(n) {
      simulate.seqs(n, locus_attrs)
    })
    names(seqs) <- rownames(locus_attrs)
    return(seqs)
  }

  # setting the set to an arbitrary value so the below is all arbitrary but
  # still deterministic.
  set.seed(0)
  # Three sets of samples across all loci
  seqs1 <- simulate.set(locus_attrs)
  seqs2 <- simulate.set(locus_attrs)
  seqs3 <- simulate.set(locus_attrs)
  set.seed(NULL)
  # 3 samples, then loci
  seqs <- list("1" = seqs1, "2" = seqs2, "3" = seqs3)

  # TODO support replicates
  write_seqs <- function(seq_sets,
                         outdir,
                         fmt="%s-%s.fasta") {
    if (! dir.exists(outdir))
      dir.create(outdir, recursive = TRUE)
    for(sn in names(seq_sets)) {
      for (ln in names(seq_sets[[sn]])) {
        fp <- file.path(outdir, sprintf(fmt, sn, ln))
        n <- names(seq_sets[[sn]][[ln]])
        if(is.null(n))
          n <- seq_along(seq_sets[[sn]][[ln]])
        dnar::write.fa(names = n,
                       dna = seq_sets[[sn]][[ln]],
                       fileName = fp)
      }
    }
  }

  prepare_for_summary <- function() {
    data.dir <- tempfile()
    write_seqs(seqs, data.dir)
    dataset <- prepare_dataset(data.dir, "()(\\d+)-([A-Za-z0-9]+).fasta")
    results <- analyze_dataset(dataset, locus_attrs,
                               summary_args = list(fraction.min = 0.05,
                                                   counts.min = 500),
                               nrepeats = 3, ncores = 1)
    lapply(dataset$Filename, file.remove)
    file.remove(data.dir)
    return(list(dataset = dataset, results = results))
  }

  results_summary_data <- prepare_for_summary()

})
