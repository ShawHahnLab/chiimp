# simulated data for testing ----------------------------------------------


# Note that having an object stored directly in the package like this (and
# forcing it to be the last file loaded when building) isn't ideal since it gets
# created and stored at build time even though the code is mixed in with the
# regular R functions.  A better way would be to explicitly build the test_data
# list and store it in data/ as Hadley describes:
# http://r-pkgs.had.co.nz/data.html

#' Helper Data for Examples
#'
#' This list is a bundle of shared data and functions for CHIIMP examples.
#' @export
test_data <- within(list(), {
  # This is a particularly awkward approach now that in the development branch
  # for version 3.6.0 the random number generator has changed its behavior.
  # The below is a stopgap measure but this should really be reorganized to not
  # need to generate the test data at build-time.
  # From ?RNGkind:
  # > sample.kind can be "Rounding" or "Rejection", or partial matches to these.
  # > The former was the default in versions prior to 3.6.0: it made sample
  # > noticeably non-uniform on large populations, and should only be used for
  # > reproduction of old results. See PR#17494 for a discussion."
  # Older R doesn't have have a third argument to RNGkind, so, only run this
  # if needed.  I'll temporarily disable warnings here so that R doesn't warn
  # about the Rounding option's behavior.
  rng_orig <- RNGkind()
  if (length(rng_orig) > 2) {
    warn_orig <- options()$warn
    options(warn = -1)
    RNGkind("Mersenne-Twister", "Inversion", "Rounding")
    options(warn = warn_orig)
    rm(warn_orig)
  }
  # Careful!  When running via a package check we might be in temporary
  # installed copy in /tmp or elsewhere, and probably won't have the "inst"
  # directory anymore.  Alternatively when running with devtools::test() we
  # will.
  f.locus_attrs <- unique(system.file(c("inst/example_locus_attrs.csv",
                                        "example_locus_attrs.csv"),
                                      package = methods::getPackageName()))
  txt.locus_attrs <- readChar(f.locus_attrs,
                              nchars = file.info(f.locus_attrs)$size)
  locus_attrs <- read.table(f.locus_attrs,
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            sep = ",")
  rm(f.locus_attrs)
  rownames(locus_attrs) <- locus_attrs$Locus

  sample.data.cols <- c("Seq", "Count", "Length", "MatchingLocus", "MotifMatch",
                        "LengthMatch", "Ambiguous", "Stutter", "Artifact",
                        "FractionOfTotal", "FractionOfLocus")
  sample.summary.cols <- c("Allele1Seq", "Allele1Count",
                           "Allele1Length", "Allele2Seq",
                           "Allele2Count", "Allele2Length",
                           "Homozygous", "Ambiguous", "Stutter", "Artifact",
                           "CountTotal", "CountLocus", "ProminentSeqs")

  make.seq_junk <- function(N) {
    nucleotides <- c("A", "T", "C", "G")
    vapply(runif(N, min = 1, max = 20), function(L)
      paste0(sample(nucleotides, L, replace = TRUE), collapse = ""),
      "character")
  }

  simulate.seqs <- function(locus_name, locus_attrs, homozygous=NULL, N=5000,
                            off_target_ratio=10, cross_contam_ratio=100) {
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
    if (off_target_ratio > 0) {
      idx <- seq(1, length(seqs), off_target_ratio)
      seqs[idx] <- paste0(attrs$Primer, make.seq_junk(10), attrs$ReversePrimer)
    }
    ## Mix in other loci
    if (cross_contam_ratio > 0) {
      others <- locus_attrs[-match(locus_name, rownames(locus_attrs)), ]
      for (i in 1:nrow(others)) {
        idx <- seq(i, length(seqs), cross_contam_ratio * nrow(others))
        seqs[idx] <- simulate.seqs(locus_name = others$Locus[i],
                                   locus_attrs = locus_attrs,
                                   N = length(idx),
                                   off_target_ratio = 0,
                                   cross_contam_ratio = 0)
      }
    }
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
    for (sn in names(seq_sets)) {
      for (ln in names(seq_sets[[sn]])) {
        fp <- file.path(outdir, sprintf(fmt, sn, ln))
        n <- names(seq_sets[[sn]][[ln]])
        if (is.null(n))
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
    results <- analyze_dataset(dataset, locus_attrs, nrepeats = 3, ncores = 1,
                               analysis_opts = list(fraction.min = 0.05),
                               summary_opts = list(counts.min = 500))
    unlink(data.dir, recursive = TRUE)
    return(list(dataset = dataset, results = results))
  }

  results_summary_data <- prepare_for_summary()

  kg1 <- subset(results_summary_data$results$summary,
         Sample == 1)[, c("Locus", "Allele1Seq", "Allele2Seq")]
  kg1 <- cbind(Name = "ID002", kg1)
  kg2 <- subset(results_summary_data$results$summary,
                Sample == 2)[, c("Locus", "Allele1Seq", "Allele2Seq")]
  kg2 <- cbind(Name = "ID001", kg2)
  genotypes_known <- rbind(kg2, kg1)

  # reset the RNG behavior
  do.call(RNGkind, as.list(rng_orig))
  rm(rng_orig)
})
