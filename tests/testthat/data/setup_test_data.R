# Use this to generate the files under tests/testthat/data/
# Run from top level of repo with the package loaded

# common ------------------------------------------------------------------

DIRPATH = "tests/testthat/data/"

# from tests/testthat/helper.R
within_tmpdir <- function(expr) {
  here <- getwd()
  data.dir <- tempfile()
  dir.create(data.dir)
  setwd(data.dir)
  tryCatch(eval(expr), finally = {
    unlink(x = data.dir, recursive = TRUE)
    setwd(here)
  })
}

# given obj or parent$obj in a setup function below, save it to
# <dirpath>/obj.rds
mktestrds <- function(obj, objname=NULL, dirpath=NULL) {
  path <- mktestfile(obj, objname, dirpath, ".rds")
  message(paste("saving", path))
  saveRDS(obj, path)
}

# given obj or parent$obj in a setup function below, save it to
# <dirpath>/obj.csv
mktestcsv <- function(obj, objname=NULL, dirpath=NULL) {
  path <- mktestfile(obj, objname, dirpath, ".csv")
  message(paste("saving", path))
  write.csv(obj, path, row.names = FALSE, quote = FALSE)
}

# don't look at this, it's horrible
mktestfile <- function(obj, objname=NULL, dirpath=NULL, ext="") {
  if (is.null(dirpath)) {
    funcname <- as.character(sys.call(-2))
    dirpath <- sub("^setup_test_data_", DIRPATH, funcname)
  }
  if (is.null(objname)) {
    objname <- deparse(substitute(obj, sys.frame(-1)))
    objname <- sub(".*\\$", "", objname)
  }
  file.path(dirpath, paste0(objname, ext))
}

# This was split away from R/zz_helper_data.R.

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

make_test_data <- function() {
  within(list(), {
    locus_attrs <- load_locus_attrs("inst/example_locus_attrs.csv")
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
    rm(seqs1)
    rm(seqs2)
    rm(seqs3)
    
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
      within_tmpdir({
        write_seqs(seqs, "data")
        dataset <- prepare_dataset("data", "()(\\d+)-([A-Za-z0-9]+).fasta")
        results <- analyze_dataset(
          dataset, locus_attrs, nrepeats = 3, ncores = 1,
          analysis_opts = list(fraction.min = 0.05),
          summary_opts = list(counts.min = 500))
      })
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
    rm(kg1)
    rm(kg2)
    
    # reset the RNG behavior
    do.call(RNGkind, as.list(rng_orig))
    rm(rng_orig)
  })
}

test_data_for_setup <- make_test_data()


# io ----------------------------------------------------------------------


setup_test_data_io <- function() {
  dirpath <- file.path(DIRPATH, "io")
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # common
  file.copy(
    "inst/example_locus_attrs.csv",
    file.path(dirpath, "locus_attrs.csv"),
    overwrite = TRUE)
  locus_attrs <- load_locus_attrs(file.path(dirpath, "locus_attrs.csv"))
  mktestrds(locus_attrs)
  # load_config
  write(
    'fp_dataset: "samples.csv"\noutput:\n  fp_rds: "results.rds"',
    file.path(dirpath, "config.yml"))
  # load config (unexpected entries)
  write(
    'fp_dataset: "samples.csv"\noutput:\n  fp_rds: "results.rds"\nunrecognized: 10\ndataset_analysis:\n  name_args:\n    unknown: 5',
    file.path(dirpath, "config_unrecognized_key.yml"))
  # load_csv (unknown columns)
  write('Vec1,Vec2,Vec3\nA,B,C\nD,E,F', file.path(dirpath, "misc.csv"))
  misc <- read.csv(file.path(dirpath, "misc.csv"))
  rownames(misc) <- paste0("entry", 1:2)
  mktestrds(misc, "misc.csv")
  # load_locus_attrs (wrong col)
  locus_attrs_wrongcol <- locus_attrs
  colnames(locus_attrs_wrongcol)[
    match("LengthMin", colnames(locus_attrs_wrongcol))] <- "length_min"
  mktestcsv(locus_attrs_wrongcol)
  # load_locus_attrs (dups)
  locus_attrs_dups <- locus_attrs
  locus_attrs_dups$Locus[2] <- "A"
  rownames(locus_attrs_dups) <- c("A", "A.1", "1", "2")
  mktestcsv(locus_attrs_dups)
  mktestrds(locus_attrs_dups)
  # load_dataset
  dataset <- expand.grid(1:3, 1:5, c("1", "2", "A", "B"))
  dataset <- data.frame(
    Filename = paste0(do.call(paste, c(dataset, list(sep = "-"))), ".fasta"),
    Replicate = as.character(dataset$Var1),
    Sample = as.character(dataset$Var2),
    Locus = as.character(dataset$Var3),
    row.names = do.call(paste, c(dataset[, c(2, 1, 3)], list(sep = "-"))),
    stringsAsFactors = FALSE)
  mktestrds(dataset)
  mktestcsv(dataset)
  # load_dataset (dups)
  extras <- subset(dataset, Sample == 1 & Replicate == 1)
  extras$Filename <- gsub(".fasta", "-alt.fasta", extras$Filename)
  rownames(extras) <- paste0(rownames(extras), ".1")
  dataset_dups <- rbind(dataset, extras)
  mktestrds(dataset_dups)
  mktestcsv(dataset_dups)
  # save_seqfile_data
  mktestrds(test_data_for_setup$seqs)
}


# analyze_seqs ------------------------------------------------------------


setup_test_data_analyze_seqs <- function() {
  dirpath <- file.path(DIRPATH, "analyze_seqs")
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # common
  locus_attrs <- load_locus_attrs("inst/example_locus_attrs.csv")
  mktestrds(locus_attrs)
  seqs <- test_data_for_setup$seqs$`1`$A
  mktestrds(seqs)
  # analyze_seqs
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  mktestrds(seq_data)
  # analyze_seqs (empty)
  seq_data_empty <- analyze_seqs(c(), locus_attrs, 3)
  mktestrds(seq_data_empty)
  # analyze_seqs (stubs)
  seqs_stubs <- seqs
  seqs_stubs[1:100] <- ""
  mktestrds(seqs_stubs)
  seq_data_stubs <- analyze_seqs(seqs_stubs, locus_attrs, 3)
  mktestrds(seq_data_stubs)
  # analyze_seqs (reverse primers)
  seq_data_rev_primers <- analyze_seqs(
    seqs, locus_attrs, 3, use_reverse_primers = TRUE)
  mktestrds(seq_data_rev_primers)
  locus_attrs_mod <- locus_attrs
  # replace locus A's reverse primer with locus B's
  locus_attrs_mod$ReversePrimer[1] <- locus_attrs_mod$ReversePrimer[2]
  mktestrds(locus_attrs_mod, "locus_attrs_rev_primer_mod")
  seq_data_rev_primers_mod <- analyze_seqs(
    seqs, locus_attrs_mod, 3, use_reverse_primers = TRUE)
  mktestrds(seq_data_rev_primers_mod)
  # analyze_seqs (reverse primers with revcmp)
  locus_attrs_revcmp <- locus_attrs
  locus_attrs_revcmp$ReversePrimer <- as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(
      locus_attrs_revcmp$ReversePrimer)))
  mktestrds(locus_attrs_revcmp, "locus_attrs_rev_primer_revcmp")
  seq_data_rev_primers_revcmp <- analyze_seqs(
    seqs, locus_attrs_revcmp, 3,
    use_reverse_primers = TRUE, reverse_primer_r1 = FALSE)
  # This should actually be the same output as earlier though
  if (! identical(seq_data_rev_primers, seq_data_rev_primers_revcmp)) {
    stop("mismatch in test data")
  }
  # analyze_seqs (high counts stutter check)
  seqs_stutter <- test_data_for_setup$seqs$`1`$A
  seqs_stutter[nchar(seqs_stutter) %in% c(158, 54)] <- seqs_stutter[
    nchar(seqs_stutter) == 190][1]
  mktestrds(seqs_stutter, "seqs_stutter_filter_check")
  seq_data_stutter <- analyze_seqs(seqs_stutter, locus_attrs, 3)
  mktestrds(seq_data_stutter, "seq_data_stutter_filter_check")
  # analyze_seqs (stutter threshold)
  seqs_stutter_thresh <- test_data_for_setup$seqs$`1`$A
  seqs_stutter_thresh[
    nchar(seqs_stutter_thresh) %in% c(158, 54)] <- seqs_stutter_thresh[
      nchar(seqs_stutter_thresh) == 190][1]
  mktestrds(seqs_stutter, "seqs_stutter_threshold_check")
  seq_data_stutter_thresh <- analyze_seqs(seqs_stutter_thresh, locus_attrs, 3)
  seq_data_stutter_thresh_mod <- analyze_seqs(
    seqs_stutter_thresh, locus_attrs, 3, stutter.count.ratio_max = 1 / 2)
  mktestrds(seq_data_stutter_thresh, "seq_data_stutter_threshold_orig")
  mktestrds(seq_data_stutter_thresh_mod, "seq_data_stutter_threshold_mod")
  # analyze_seqs (artifact)
  seqs_artifact <- test_data_for_setup$seqs$`1`$A
  # Take that first stutter and make it an artifact instead
  highest <- names(sort(table(seqs_artifact), decreasing = TRUE)[1])
  stutter <- names(sort(table(seqs_artifact), decreasing = TRUE)[3])
  idx <- seqs_artifact == stutter
  seqs_artifact[idx] <- highest
  substr(seqs_artifact[idx], nchar(stutter), nchar(stutter)) <- "R"
  seq_data_artifact <- analyze_seqs(seqs_artifact, locus_attrs, 3)
  mktestrds(seqs_artifact)
  mktestrds(seq_data_artifact)
  # analyze_seqs (ambig)
  seqs_ambig <- test_data_for_setup$seqs$`1`$A
  seqs_ambig[seqs_ambig == seqs_ambig[1]] <- sub(
    "AGCCAGTC", "AGCCANTC", seqs_ambig[1])
  seq_data_ambig <- analyze_seqs(seqs_ambig, locus_attrs, 3)
  mktestrds(seqs_ambig)
  mktestrds(seq_data_ambig)
}


# analyze_sample ----------------------------------------------------------


setup_test_data_analyze_sample <- function() {
  dirpath <- file.path(DIRPATH, "analyze_sample")
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # analyze_sample
  seq_data <- with(
    test_data_for_setup, analyze_seqs(seqs$`1`$A, locus_attrs, 3))
  sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
  mktestrds(seq_data)
  mktestrds(sample_data)
}


# categorize --------------------------------------------------------------


setup_test_data_categorize <- function(
  dirpath="tests/testthat/data/categorize") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # match_known_genotypes
  # Names for the first two samples but leaving the third unidentified
  results_summary <- test_data_for_setup$results_summary_data$results$summary
  results_summary$Name <- c("ID002", "ID001")[
    as.integer(results_summary$Sample)]
  mktestrds(results_summary)
  mktestrds(test_data_for_setup$genotypes_known)
  genotypes_matched <- with(test_data_for_setup, match_known_genotypes(
    results_summary, genotypes_known))
  mktestrds(genotypes_matched, "match_known_genotypes")
  # categorize_genotype_results
  results_summary_matched <- cbind(results_summary, genotypes_matched)
  mktestrds(results_summary_matched)
  categories <- categorize_genotype_results(results_summary_matched)
  mktestrds(categories)
  # categorize_genotype_results (drop)
  # drop an allele for each of two samples
  results_summary_matched_drop <- within(results_summary_matched, {
    Allele1Seq[2] <- NA
    Allele2Seq[1] <- NA
  })
  categories_drop <- categorize_genotype_results(results_summary_matched_drop)
  mktestrds(results_summary_matched_drop)
  mktestrds(categories_drop)
  # categorize_genotype_results (blanks)
  results_summary_matched_blanks <- within(results_summary_matched, {
    Allele1Seq[c(1, 3)] <- NA
    Allele2Seq[1] <- NA
  })
  categories_blanks <- categorize_genotype_results(
    results_summary_matched_blanks)
  mktestrds(results_summary_matched_blanks)
  mktestrds(categories_blanks)
}


# summarize_sample --------------------------------------------------------


setup_test_data_summarize_sample <- function(
  dirpath="tests/testthat/data/summarize_sample") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # summarize_sample
  seq_data <- with(
    test_data_for_setup, analyze_seqs(seqs$`1`$A, locus_attrs, 3))
  sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
  sample_summary <- summarize_sample(
    sample_data, list(Locus = "A"), counts.min = 500)
  mktestrds(sample_data)
  mktestrds(sample_summary)
  # summarize_sample (empty)
  seq_data_empty <- analyze_seqs(c(), test_data_for_setup$locus_attrs, 3)
  sample_data_empty <- analyze_sample(seq_data_empty, list(Locus = "A"), 0.05)
  sample_summary_empty <- summarize_sample(
    sample_data_empty, list(Locus = "A"), counts.min = 500)
  mktestrds(sample_data_empty)
  mktestrds(sample_summary_empty)
  # summarize_sample (stubs)
  sample_data_stubs <- with(test_data_for_setup, {
    seqs <- test_data_for_setup$seqs$`1`$A
    seqs[1:100] <- ""
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    analyze_sample(seq_data, list(Locus = "A"), 0.05)
  })
  sample_summary_stubs <- summarize_sample(
    sample_data_stubs, list(Locus = "A"), counts.min = 500)
  mktestrds(sample_data_stubs)
  mktestrds(sample_summary_stubs)
  # summarize_sample (stutter)
  sample_data_stutter <- with(test_data_for_setup, {
    seqs <- test_data_for_setup$seqs$`3`$A
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    analyze_sample(seq_data, list(Locus = "A"), 0.05)
  })
  sample_summary_stutter <- summarize_sample(
    sample_data_stutter, list(Locus = "A"), counts.min = 500)
  mktestrds(sample_data_stutter)
  mktestrds(sample_summary_stutter)
  # summarize_sample (multi artifact)
  sample_data_multi_artifact <- with(test_data_for_setup, {
    seqs <- test_data_for_setup$seqs$`2`$`2`
    idx <- (which(seqs == seqs[1]))[c(TRUE, FALSE)] # every other matching index
    seqs[idx] <- gsub(".$", "N", seqs[idx]) # replace last character with N
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    analyze_sample(seq_data, list(Locus = "2"), 0.05)
  })
  sample_summary_multi_artifact <- summarize_sample(
    sample_data_multi_artifact, counts.min = 500)
  mktestrds(sample_data_multi_artifact)
  mktestrds(sample_summary_multi_artifact)
  # summarize_sample (multi stutter)
  sample_data_multi_stutter <- with(test_data_for_setup, {
    # Replace the third entry with a different stutter sequence.  Munge the
    # counts around to still total correctly.
    seq_data <- analyze_seqs(seqs$`3`$A, locus_attrs, 3)
    tot <- sum(seq_data$Count)
    seq_data[3, ] <- seq_data[2, ]
    seq_data[3, "Seq"] <- sub("TAGA", "TACA", seq_data[3, "Seq"])
    seq_data[3, "Count"] <- 410
    seq_data[3, "FractionOfTotal"] <- 410 / tot
    seq_data[3, "FractionOfLocus"] <- 410 / tot
    seq_data[4:12, "Count"] <- 10
    seq_data[4:12, "FractionOfTotal"] <- 10 / tot
    seq_data[4:12, "FractionOfLocus"] <- 10 / tot
    analyze_sample(seq_data, list(Locus = "A"), 0.05)
  })
  sample_summary_multi_stutter <- summarize_sample(
    sample_data_multi_stutter, counts.min = 500)
  mktestrds(sample_data_multi_stutter)
  mktestrds(sample_summary_multi_stutter)
  # summarize_sample (ambig)
  sample_data_ambig <- with(test_data_for_setup, {
    # Replace the second-highest-count sequence to include an ambiguous base
    # near the end.
    seqs <- seqs$`3`$A
    idx <- nchar(seqs) == 178
    seqs[idx] <- sub("AGCCAGTC$", "AGCCNAGTC", seqs[idx])
    seq_data <- analyze_seqs(seqs, locus_attrs, 3)
    analyze_sample(seq_data, list(Locus = "A"), 0.05)
  })
  sample_summary_ambig <- summarize_sample(sample_data_ambig, counts.min = 500)
  mktestrds(sample_data_ambig)
  mktestrds(sample_summary_ambig)
  # summarize_sample (low)
  sample_data_low <- with(test_data_for_setup, {
    # Replace the second-highest-count sequence to include an ambiguous base
    # near the end.
    seq_data <- analyze_seqs(seqs$`1`$A, locus_attrs, 3)
    seq_data$Count <- round(seq_data$Count / 100)
    analyze_sample(seq_data, list(Locus = "A"), 0.05)
  })
  sample_summary_low <- summarize_sample(sample_data_low, counts.min = 500)
  mktestrds(sample_data_low)
  mktestrds(sample_summary_low)
  # summarize_sample (prominent seqs)
  sample_datas_b <- lapply(1:3, function(samp) {
    with(test_data_for_setup, {
    seq_data <- analyze_seqs(seqs[[samp]]$B, locus_attrs, 3)
    sample_data <- analyze_sample(seq_data, list(Locus = "B"), 0.05)
    })
  })
  sample_summaries_b <- lapply(
    sample_datas_b, summarize_sample, list(Locus = "B"), counts.min = 500)
  mktestrds(sample_datas_b)
  mktestrds(sample_summaries_b)

}


# histogram ---------------------------------------------------------------


setup_test_data_histogram <- function(
  dirpath="tests/testthat/data/histogram") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # histogram
  results <- test_data_for_setup$results_summary_data$results
  seq_data <- results$files[["data/3-2.fasta"]]
  sample_data <- results$samples[["3-2"]]
  png(fp_devnull)
  output <- histogram(
    seq_data = seq_data, sample_data = sample_data, cutoff_fraction = 0.05)
  dev.off()
  mktestrds(seq_data)
  mktestrds(sample_data)
  mktestrds(output)
  # histogram (seq_data only)
  png(fp_devnull)
  output_seq_data_only <- histogram(seq_data = seq_data)
  dev.off()
  mktestrds(output_seq_data_only)
  # histogram (empty sample data)
  png(fp_devnull)
  output_empty_sample_data <- histogram(
    seq_data = seq_data, sample_data = sample_data[0, ], cutoff_fraction = 0.05)
  dev.off()
  mktestrds(output_empty_sample_data)
}


# analyze_dataset ---------------------------------------------------------


setup_test_data_analyze_dataset <- function(
  dirpath="tests/testthat/data/analyze_dataset") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # analyze_dataset
  mktestrds(test_data_for_setup$seqs)
  mktestrds(test_data_for_setup$results_summary_data$dataset)
  mktestrds(test_data_for_setup$locus_attrs)
  mktestrds(test_data_for_setup$results_summary_data$results)
  # analyze_dataset (names known alleles)
  known_alleles <- data.frame(
    Locus = c("1", "1", "A"),
    Seq = c(paste0(
      "ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTA",
      "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
      "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
      "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
      "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTGTGGCTCA",
      "AAAGCTGAAT"), paste0(
      "ACAGTCAAGAATAACTGCCCTATCTATCTATCTATCTATCTATCTATCTATCTA",
      "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
      "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
      "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC",
      "TATCTATCTATCTATCTATCTATCCTGTGGCTCAAAAGCTGAAT"), paste0(
      "TATCACTGGTGTTAGTCCTCTGTAGATAGATAGATAGATAGATAGATAGATAGA",
      "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA",
      "GATAGATAGATAGATAGATAGATAGATAGATAGACACAGTTGTGTGAGCCAGTC")),
    Name = c("280-a", "260-X", "different_name_format"))
  within_tmpdir({
    results_known_alleles <- with(test_data_for_setup, {
      write_seqs(seqs, "data")
      dataset <- prepare_dataset("data", "()(\\d+)-([A-Za-z0-9]+).fasta")
      analyze_dataset(
        dataset, locus_attrs, nrepeats = 3, ncores = 1,
        analysis_opts = list(fraction.min = 0.05),
        summary_opts = list(counts.min = 500), known_alleles = known_alleles)
    })
  })
  mktestrds(known_alleles)
  mktestrds(results_known_alleles)
}



# report ------------------------------------------------------------------


setup_test_data_report <- function(
  dirpath="tests/testthat/data/report") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  results_summary <- test_data_for_setup$results_summary_data$results$summary
  mktestrds(results_summary)
  # plot_alignment
  alignments <- align_alleles(results_summary)
  mktestrds(alignments)
  # plot_dist_mat
  dist_mat <- make_dist_mat(results_summary)
  mktestrds(dist_mat)
  # tabulate_allele_names
  allele_summary <- tabulate_allele_names(results_summary)
  mktestrds(allele_summary)
  # plot_cts_per_locus
  results <- test_data_for_setup$results_summary_data$results
  results <- summarize_dataset(results)
  mktestrds(results)
}


# all ---------------------------------------------------------------------


setup_test_data_io()
setup_test_data_analyze_seqs()
setup_test_data_analyze_sample()
setup_test_data_categorize()
setup_test_data_summarize_sample()
setup_test_data_histogram()
setup_test_data_analyze_dataset()
setup_test_data_report()