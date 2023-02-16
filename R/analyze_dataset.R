# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze_seqs and
# summarize_sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a summary data frame showing
#' each sample's summary row-by-row, a list of processed input files, and a list
#' of processed samples.  The entries in the processed-samples list and the rows
#' in the summary data frame will be sorted according to the ordering of loci in
#' `locus_attrs` and by the sample attributes.  Processed files are stored
#' separately (as there may be multiple samples per file) and named by input
#' file path.  An error is thrown if any locus entries in the given dataset are
#' not found in the locus attributes data frame.
#'
#' @param dataset data frame of sample details as produced by [prepare_dataset].
#' @param locus_attrs data frame of locus attributes as produced by
#'   [load_locus_attrs].
#' @param analysis_opts list of supplemental arguments to `analysis_function`.
#' @param summary_opts list of supplemental arguments to `summary_function`.
#' @param analysis_function function to use when analyzing each sample's data
#'   frame into the filtered version  Defaults to [analyze_sample].
#' @param summary_function function to use when summarizing each sample's full
#'   details into the standard attributes.  Defaults to [summarize_sample].
#' @param ncores integer number of CPU cores to use in parallel for sample
#'   analysis.  Defaults to one less than half the number of detected cores with
#'   a minimum of 1.  If 1, the function will run without using the
#'   [parallel] package.
#' @param known_alleles data frame of custom allele names as defined for
#'   [load_allele_names].  if NULL only the names automatically
#'   generated for the dataset summary will be used.
#'
#' @return list of results, with `summary` set to the single summary data frame,
#'   `files` the processed sequence files, and `samples` the per-sample data
#'   frames.
#'
#' @export
#' @md
analyze_dataset <- function(
    dataset,
    locus_attrs,
    analysis_opts,
    summary_opts,
    analysis_function = analyze_sample,
    summary_function = summarize_sample,
    ncores = cfg("ncores"),
    known_alleles = NULL) {
  if (! all(dataset$Locus %in% locus_attrs$Locus)) {
    rogue_loci <- unique(dataset$Locus[! dataset$Locus %in% locus_attrs$Locus])
    msg <- paste("ERROR: Locus names in dataset not in attributes table:",
                 paste(rogue_loci, collapse = ", "))
    stop(msg)
  }
  analyze_file <- function(fp, locus_attrs) {
    seqs <- load_seqs(fp)
    analyze_seqs(seqs, locus_attrs)
  }
  analyze_entry <- function(
      entry, analysis_opts, summary_opts, analysis_function, summary_function,
      analyzed_files) {
    # Get all data from the relevant file
    seq_data <- analyzed_files[[entry["Filename"]]]
    # Process into single-sample data frame
    analysis_args <- c(
      list(seq_data = seq_data, sample_attrs = entry), analysis_opts)
    sample_data <- do.call(analysis_function, analysis_args)
    # Process into single-sample summary list
    summary_args <- c(
      list(sample_data = sample_data, sample_attrs = entry), summary_opts)
    sample.summary <- do.call(summary_function, summary_args)
    # Return the processed per-sample data
    return(list(summary = sample.summary, data = sample_data))
  }
  if (ncores > 1) {
    # Set up the cluster and export required names (those objects used in
    # analyze_entry that are expected from the environment and not passed as
    # arguments).
    cluster_names <- c("locus_attrs")
    cluster <- parallel::makeCluster(ncores)
    # https://stackoverflow.com/a/12232695/6073858
    parallel::clusterEvalQ(cluster, library(dnar))
    # Load the currently-used version of chiimp if running from a source
    # version.
    # Without this, parallel runs will use an installed version of chiimp, if
    # present, even if the first process started from a devtools source
    # directory.  This would give unexpected behavior between different
    # installed/source versions.  There's probably a better way to handle this
    # situation but this works for now.
    pkg <- methods::getPackageName()
    if (basename(system.file(package = pkg)) == "inst") {
      dp <- dirname(system.file(package = pkg))
      parallel::clusterExport(cluster, "dp", envir = environment())
      parallel::clusterEvalQ(cluster, devtools::load_all(dp))
    }
    parallel::clusterExport(
      cl = cluster, varlist = cluster_names, envir = environment())
    tryCatch({
      # Load, analyze, and summarize each sample across the cluster.  Each row
      # in the dataset data frame will be given as the entry argument to
      # analyze_entry.
      fps <- unique(dataset$Filename)
      analyzed_files <- parallel::parLapply(
        cluster, fps,
        analyze_file, locus_attrs = locus_attrs)
      names(analyzed_files) <- fps
      raw_results <- parallel::parApply(
        cluster, dataset, 1, analyze_entry,
        analysis_opts = analysis_opts,
        summary_opts = summary_opts,
        analysis_function = analysis_function,
        summary_function = summary_function,
        analyzed_files = analyzed_files)
    },
    finally = {
      parallel::stopCluster(cluster)
    })
  } else {
    fps <- unique(dataset$Filename)
    analyzed_files <- lapply(
      fps,
      analyze_file, locus_attrs = locus_attrs)
    names(analyzed_files) <- fps
    raw_results <- apply(
      dataset, 1, analyze_entry,
      analysis_opts = analysis_opts,
      summary_opts = summary_opts,
      analysis_function = analysis_function,
      summary_function = summary_function,
      analyzed_files = analyzed_files)
  }

  # Check if any of the raw data files had no reads to start with.
  empties <- sum(sapply(analyzed_files, nrow) == 0)
  if (empties) {
    logmsg(paste("WARNING: Zero reads for", empties, "of",
                 length(analyzed_files), "data files"))
  }

  # Recombined results into a summary data frame and a list of full sample data.
  results <- tidy_analyzed_dataset(dataset, raw_results)
  results$files <- analyzed_files
  # Add allele name columns to all data frames for any allele in the given
  # known_alleles data frame or called in the current genotypes.
  results <- name_known_sequences(results, known_alleles)
  # Reorder entries to match locus_attrs.
  results <- sort_results(results, locus_attrs)
  results
}

#' Tidy raw analyzed dataset results
#'
#' Rearrange pairs of sample summary / sample data objects into a single
#' summary cross-sample data frame and a list of detailed per-sample data
#' frames.
#'
#' @param dataset data frame of sample details as produced by [prepare_dataset].
#' @param raw_results list of pairs of sample summary and sample data data
#'   frames (from [summarize_sample]) and [analyze_seqs]).
#'
#' @return list of results, with `summary` set to the single summary data frame
#'   and `data` the per-sample data frames.
#' @md
tidy_analyzed_dataset <- function(dataset, raw_results) {
  summaries <- lapply(raw_results, function(s) {
    data.frame(s[[1]], stringsAsFactors = FALSE)
  })
  samples <- lapply(raw_results, `[[`, 2)
  summary <- do.call(rbind, summaries)
  rownames(summary) <- rownames(dataset)
  names(samples)    <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, samples = samples))
}

#' Name known allele sequences
#'
#' For the given results list (pair of summary data frame and list of per-sample
#' data frames as produced by [tidy_analyzed_dataset]), add columns
#' to all data frames defining names for recognized sequences.  For the summary
#' data frame this will be `Allele1Name` and `Allele2Name`.  For each
#' sample data frame this will be `SeqName`, defined for any sequences
#' represented in the summary or in a given known alleles set.
#'
#' @param results results list as produced by `tidy_analyzed_dataset`.
#' @param known_alleles data frame of custom allele names as defined for
#'   `load_allele_names`.  if NULL only the names automatically generated for
#'   the summary will be used.
#' @param ... additional arguments to `make_allele_names`.
#'
#' @return list of results, with `summary` set to the single summary data frame
#'   and `data` the per-sample data frames.  A `SeqName` column in sample data
#'   frames and `Allele1Name` and `Allele2Name` columns in the summary data
#'   frame will associate any sequence matching a known allele (for either the
#'   given table or the current dataset) with a text name.
#' @md
name_known_sequences <- function(results, known_alleles, ...) {
  # Name all of the called alleles across samples
  results$summary <- name_alleles_in_table(
    results$summary, known_alleles, ...)

  # Create table of allele names for current dataset
  a1 <- results$summary[, c("Locus", "Allele1Seq", "Allele1Name")]
  a2 <- results$summary[, c("Locus", "Allele2Seq", "Allele2Name")]
  colnames(a1) <- c("Locus", "Seq", "Name")
  colnames(a2) <- c("Locus", "Seq", "Name")
  a_combo <- rbind(a1, a2)
  a_combo <- unique(a_combo[! is.na(a_combo$Seq), ])

  # Merge into given known alleles table (if present)
  known_alleles <- if (is.null(known_alleles)) {
    a_combo
  } else {
    known_alleles <- known_alleles[, c("Locus", "Seq", "Name")]
    known_alleles$Name <- as.character(known_alleles$Name)
    unique(rbind(known_alleles, a_combo))
  }

  # Name recognized sequences in each sample data frame
  results$samples <- lapply(results$samples, function(d) {
    idx <- match(d$Seq, known_alleles$Seq)
    d$SeqName <- known_alleles$Name[idx]
    d
  })

  return(results)
}

#' Sort entries in results data frames
#'
#' Rearrange rows in the summary data frame and corresponding entries in the
#' per-sample data frames list, matching locus order given in locus attributes.
#' The Locus column of the summary data frame will be set to a factor to
#' preserve the defined order.  Only levels remaining in that set are kept.
#'
#' @param results results list as produced by `tidy_analyzed_dataset`.
#' @param locus_attrs data frame of locus attributes as produced by
#'   `load_locus_attrs`.
#'
#' @return list of results, with `summary` set to the single summary data frame
#'   and `data` the per-sample data frames.  `summary$Locus` is coerced to a
#'   factor with levels ordered according to their appearance in
#'   `locus_attrs$Locus`.  Order of rows in `summary` and entries in `data` are
#'   updated accordingly.
#' @md
sort_results <- function(results, locus_attrs) {
  results$summary$Locus <- factor(results$summary$Locus,
                                  levels = locus_attrs$Locus)
  results$summary$Locus <- droplevels(results$summary$Locus)
  ord <- order_entries(results$summary)
  results$summary <- results$summary[ord, ]
  results$samples <- results$samples[ord]
  results
}
