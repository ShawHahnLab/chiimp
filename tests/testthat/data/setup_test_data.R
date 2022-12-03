# Use this to generate the files under tests/testthat/data/
# Run from top level of repo with the package loaded

# common ------------------------------------------------------------------


# given obj or parent$obj in a setup function below, save it to
# <dirpath>/obj.rds
mktestrds <- function(obj, objname=NULL) {
  dirpath <- get("dirpath", env=sys.frame(-1))
  if (is.null(objname)) {
    objname <- sub(".*\\$", "", deparse(substitute(obj)))
  }
  path <- file.path(dirpath, paste0(objname, ".rds"))
  message(paste("saving", path))
  saveRDS(obj, path)
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
    rm(kg1)
    rm(kg2)
    
    # reset the RNG behavior
    do.call(RNGkind, as.list(rng_orig))
    rm(rng_orig)
  })
}

test_data_for_setup <- make_test_data()


# io ----------------------------------------------------------------------


setup_test_data_io <- function(dirpath="tests/testthat/data/io") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # common
  file.copy(
    "inst/example_locus_attrs.csv",
    file.path(dirpath, "locus_attrs.csv"),
    overwrite = TRUE)
  locus_attrs <- load_locus_attrs(file.path(dirpath, "locus_attrs.csv"))
  saveRDS(locus_attrs, file.path(dirpath, "locus_attrs.rds"))
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
  saveRDS(misc, file.path(dirpath, "misc.csv.rds"))
  # load_locus_attrs (wrong col)
  locus_attrs_wrongcol <- locus_attrs
  colnames(locus_attrs_wrongcol)[
    match("LengthMin", colnames(locus_attrs_wrongcol))] <- "length_min"
  write.csv(
    locus_attrs_wrongcol, file.path(dirpath, "locus_attrs_wrongcol.csv"),
    row.names = FALSE, quote = FALSE)
  # load_locus_attrs (dups)
  locus_attrs_dups <- locus_attrs
  locus_attrs_dups$Locus[2] <- "A"
  rownames(locus_attrs_dups) <- c("A", "A.1", "1", "2")
  write.csv(
    locus_attrs_dups, file.path(dirpath, "locus_attrs_dups.csv"),
    row.names = FALSE, quote = FALSE)
  saveRDS(locus_attrs_dups, file.path(dirpath, "locus_attrs_dups.rds"))
  # load_dataset
  dataset <- expand.grid(1:3, 1:5, c("1", "2", "A", "B"))
  dataset <- data.frame(
    Filename = paste0(do.call(paste, c(dataset, list(sep = "-"))), ".fasta"),
    Replicate = as.character(dataset$Var1),
    Sample = as.character(dataset$Var2),
    Locus = as.character(dataset$Var3),
    row.names = do.call(paste, c(dataset[, c(2, 1, 3)], list(sep = "-"))),
    stringsAsFactors = FALSE)
  saveRDS(dataset, file.path(dirpath, "dataset.rds"))
  write.csv(
    dataset, file.path(dirpath, "dataset.csv"),
    row.names = FALSE, quote = FALSE)
  # load_dataset (dups)
  extras <- subset(dataset, Sample == 1 & Replicate == 1)
  extras$Filename <- gsub(".fasta", "-alt.fasta", extras$Filename)
  rownames(extras) <- paste0(rownames(extras), ".1")
  dataset_dups <- rbind(dataset, extras)
  saveRDS(dataset_dups, file.path(dirpath, "dataset_dups.rds"))
  write.csv(
    dataset_dups, file.path(dirpath, "dataset_dups.csv"),
    row.names = FALSE, quote = FALSE)
  # save_seqfile_data
  saveRDS(test_data_for_setup$seqs, file.path(dirpath, "seqs.rds"))
}


# analyze_seqs ------------------------------------------------------------


setup_test_data_analyze_seqs <- function(
  dirpath="tests/testthat/data/analyze_seqs") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # common
  locus_attrs <- load_locus_attrs("inst/example_locus_attrs.csv")
  saveRDS(locus_attrs, file.path(dirpath, "locus_attrs.rds"))
  seqs <- test_data_for_setup$seqs$`1`$A
  saveRDS(seqs, file.path(dirpath, "seqs.rds"))
  # analyze_seqs
  seq_data <- analyze_seqs(seqs, locus_attrs, 3)
  saveRDS(seq_data, file.path(dirpath, "seq_data.rds"))
  # analyze_seqs (empty)
  seq_data_empty <- analyze_seqs(c(), locus_attrs, 3)
  saveRDS(seq_data_empty, file.path(dirpath, "seq_data_empty.rds"))
  # analyze_seqs (stubs)
  seqs_stubs <- seqs
  seqs_stubs[1:100] <- ""
  saveRDS(seqs_stubs, file.path(dirpath, "seqs_stubs.rds"))
  seq_data_stubs <- analyze_seqs(seqs_stubs, locus_attrs, 3)
  saveRDS(seq_data_stubs, file.path(dirpath, "seq_data_stubs.rds"))
  # analyze_seqs (reverse primers)
  seq_data_rev_primers <- analyze_seqs(
    seqs, locus_attrs, 3, use_reverse_primers = TRUE)
  saveRDS(seq_data_rev_primers, file.path(dirpath, "seq_data_rev_primers.rds"))
  locus_attrs_mod <- locus_attrs
  # replace locus A's reverse primer with locus B's
  locus_attrs_mod$ReversePrimer[1] <- locus_attrs_mod$ReversePrimer[2]
  saveRDS(
    locus_attrs_mod, file.path(dirpath, "locus_attrs_rev_primer_mod.rds"))
  seq_data_rev_primers_mod <- analyze_seqs(
    seqs, locus_attrs_mod, 3, use_reverse_primers = TRUE)
  saveRDS(
    seq_data_rev_primers_mod,
    file.path(dirpath, "seq_data_rev_primers_mod.rds"))
  # analyze_seqs (reverse primers with revcmp)
  locus_attrs_revcmp <- locus_attrs
  locus_attrs_revcmp$ReversePrimer <- as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(
      locus_attrs_revcmp$ReversePrimer)))
  saveRDS(
    locus_attrs_revcmp,
    file.path(dirpath, "locus_attrs_rev_primer_revcmp.rds"))
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
  saveRDS(seqs_stutter, file.path(dirpath, "seqs_stutter_filter_check.rds"))
  seq_data_stutter <- analyze_seqs(seqs_stutter, locus_attrs, 3)
  saveRDS(
    seq_data_stutter,
    file.path(dirpath, "seq_data_stutter_filter_check.rds"))
  # analyze_seqs (stutter threshold)
  seqs_stutter_thresh <- test_data_for_setup$seqs$`1`$A
  seqs_stutter_thresh[
    nchar(seqs_stutter_thresh) %in% c(158, 54)] <- seqs_stutter_thresh[
      nchar(seqs_stutter_thresh) == 190][1]
  saveRDS(seqs_stutter, file.path(dirpath, "seqs_stutter_threshold_check.rds"))
  seq_data_stutter_thresh <- analyze_seqs(seqs_stutter_thresh, locus_attrs, 3)
  seq_data_stutter_thresh_mod <- analyze_seqs(
    seqs_stutter_thresh, locus_attrs, 3, stutter.count.ratio_max = 1 / 2)
  saveRDS(
    seq_data_stutter_thresh,
    file.path(dirpath, "seq_data_stutter_threshold_orig.rds"))
  saveRDS(
    seq_data_stutter_thresh_mod,
    file.path(dirpath, "seq_data_stutter_threshold_mod.rds"))
  # analyze_seqs (artifact)
  seqs_artifact <- test_data_for_setup$seqs$`1`$A
  # Take that first stutter and make it an artifact instead
  highest <- names(sort(table(seqs_artifact), decreasing = TRUE)[1])
  stutter <- names(sort(table(seqs_artifact), decreasing = TRUE)[3])
  idx <- seqs_artifact == stutter
  seqs_artifact[idx] <- highest
  substr(seqs_artifact[idx], nchar(stutter), nchar(stutter)) <- "R"
  seq_data_artifact <- analyze_seqs(seqs_artifact, locus_attrs, 3)
  saveRDS(seqs_artifact, file.path(dirpath, "seqs_artifact.rds"))
  saveRDS(seq_data_artifact, file.path(dirpath, "seq_data_artifact.rds"))
  # analyze_seqs (ambig)
  seqs_ambig <- test_data_for_setup$seqs$`1`$A
  seqs_ambig[seqs_ambig == seqs_ambig[1]] <- sub(
    "AGCCAGTC", "AGCCANTC", seqs_ambig[1])
  seq_data_ambig <- analyze_seqs(seqs_ambig, locus_attrs, 3)
  saveRDS(seqs_ambig, file.path(dirpath, "seqs_ambig.rds"))
  saveRDS(seq_data_ambig, file.path(dirpath, "seq_data_ambig.rds"))
}


# analyze_sample ----------------------------------------------------------


setup_test_data_analyze_sample <- function(
  dirpath="tests/testthat/data/analyze_sample") {
  if (! dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }
  # analyze_sample
  seq_data <- with(
    test_data_for_setup, analyze_seqs(seqs$`1`$A, locus_attrs, 3))
  sample_data <- analyze_sample(seq_data, list(Locus = "A"), 0.05)
  saveRDS(seq_data, file.path(dirpath, "seq_data.rds"))
  saveRDS(sample_data, file.path(dirpath, "sample_data.rds"))
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


# all ---------------------------------------------------------------------


setup_test_data_io()
setup_test_data_analyze_seqs()
setup_test_data_analyze_sample()
setup_test_data_categorize()
