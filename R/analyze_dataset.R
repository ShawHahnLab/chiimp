# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze_seqs and
# summarize_sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a summary data frame showing
#' each sample's summary row-by-row, a list of processed input files, and a list
#' of processed samples.  The entries in the processed-samples list and the rows
#' in the summary data frame will be sorted according to the ordering of loci in
#' \code{locus_attrs} and by the sample attributes.  Processed files are stored
#' separately (as there may be multiple samples per file) and named by input
#' file path.  An error is thrown if any locus entries in the given dataset are
#' not found in the locus attributes data frame.
#'
#' @param dataset data frame of sample details as produced by
#'   \code{\link{prepare_dataset}}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{\link{load_locus_attrs}}.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match (see \code{\link{analyze_seqs}}).
#' @param ncores integer number of CPU cores to use in parallel for sample
#'   analysis.  Defaults to one less than half the number of detected cores with
#'   a minimum of 1.  If 1, the function will run without using the
#'   \code{parallel} package.
#' @param analysis_opts list of supplemental arguments to
#'   \code{analysis_function}.
#' @param summary_opts list of supplemental arguments to
#'   \code{summary_function}.
#' @param analysis_function function to use when analyzing each sample's data
#'   frame into the filtered version  Defaults to \code{\link{analyze_sample}}.
#' @param summary_function function to use when summarizing each sample's full
#'   details into the standard attributes.  Defaults to
#'   \code{\link{summarize_sample}}.
#' @param known_alleles data frame of custom allele names as defined for
#'   \code{\link{load_allele_names}}.  if NULL only the names automatically
#'   generated for the dataset summary will be used.
#' @param name_args list of additional arguments to
#'   \code{\link{make_allele_name}}.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame, \code{files} the processed sequence files, and \code{samples} the
#'   per-sample data frames.
#'
#' @export
analyze_dataset <- function(dataset,
                            locus_attrs,
                            nrepeats,
                            ncores = 0,
                            analysis_opts,
                            summary_opts,
                            analysis_function=analyze_sample,
                            summary_function=summarize_sample,
                            known_alleles=NULL,
                            name_args=list()) {
  if (! all(dataset$Locus %in% locus_attrs$Locus)) {
    rogue_loci <- unique(dataset$Locus[! dataset$Locus %in% locus_attrs$Locus])
    msg <- paste("ERROR: Locus names in dataset not in attributes table:",
                 paste(rogue_loci, collapse = ", "))
    stop(msg)
  }
  if (ncores == 0) {
    ncores <- max(1, as.integer(parallel::detectCores() / 2) - 1)
  }
  analyze.file <- function(fp, locus_attrs, nrepeats) {
    seqs <- load_seqs(fp)
    analyze_seqs(seqs, locus_attrs, nrepeats)
  }
  analyze.entry <- function(entry, analysis_opts, summary_opts,
                            analysis_function, summary_function,
                            analyzed_files) {
    # Get all data from the relevant file
    seq_data <- analyzed_files[[entry["Filename"]]]
    # Process into single-sample data frame
    analysis_args <- c(list(seq_data = seq_data, sample.attrs = entry),
                       analysis_opts)
    sample_data <- do.call(analysis_function, analysis_args)
    # Process into single-sample summary list
    summary_args <- c(list(sample_data = sample_data, sample.attrs = entry),
                      summary_opts)
    sample.summary <- do.call(summary_function, summary_args)
    # Return the processed per-sample data
    return(list(summary = sample.summary, data = sample_data))
  }
  if (ncores > 1) {
    # Set up the cluster and export required names (those objects used in
    # analyze.entry that are expected from the environment and not passed as
    # arguments).
    cluster_names <- c("locus_attrs", "nrepeats")
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
    parallel::clusterExport(cl = cluster,
                            varlist = cluster_names,
                            envir = environment())
    tryCatch({
      # Load, analyze, and summarize each sample across the cluster.  Each row
      # in the dataset data frame will be given as the entry argument to
      # analyze.entry.
      fps <- unique(dataset$Filename)
      analyzed_files <- parallel::parLapply(cluster, fps, analyze.file,
                                            locus_attrs = locus_attrs,
                                            nrepeats = nrepeats)
      names(analyzed_files) <- fps
      raw.results <- parallel::parApply(cluster, dataset, 1, analyze.entry,
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
    analyzed_files <- lapply(fps, analyze.file,
                             locus_attrs = locus_attrs,
                             nrepeats = nrepeats)
    names(analyzed_files) <- fps
    raw.results <- apply(dataset, 1, analyze.entry,
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
  results <- tidy_analyzed_dataset(dataset, raw.results)
  results$files <- analyzed_files
  # Add allele name columns to all data frames for any allele in the given
  # known_alleles data frame or called in the current genotypes.
  results <- name_known_sequences(results, known_alleles, name_args)
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
#' @param dataset data frame of sample details as produced by
#'   \code{\link{prepare_dataset}}.
#' @param raw.results list of pairs of sample summary and sample data data
#'   frames (from \code{\link{summarize_sample}} and
#'   \code{\link{analyze_seqs}}).
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.
tidy_analyzed_dataset <- function(dataset, raw.results) {
  summaries <- lapply(raw.results, function(s) {
    data.frame(s[[1]], stringsAsFactors = FALSE)
  })
  samples <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind, summaries)
  rownames(summary) <- rownames(dataset)
  names(samples)    <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, samples = samples))
}

#' Name known allele sequences
#'
#' For the given results list (pair of summary data frame and list of per-sample
#' data frames as produced by \code{\link{tidy_analyzed_dataset}}), add columns
#' to all data frames defining names for recognized sequences.  For the summary
#' data frame this will be \code{Allele1Name} and \code{Allele2Name}.  For each
#' sample data frame this will be \code{SeqName}, defined for any sequences
#' represented in the summary or in a given known alleles set.
#'
#' @param results results list as produced by
#'   \code{\link{tidy_analyzed_dataset}}.
#' @param known_alleles data frame of custom allele names as defined for
#'   \code{\link{load_allele_names}}.  if NULL only the names automatically
#'   generated for the summary will be used.
#' @param name_args list of additional arguments to
#'   \code{\link{make_allele_name}}.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.  A \code{SeqName} column
#'   in sample data frames and \code{Allele1Name} and \code{Allele2Name} columns
#'   in the summary data frame will associate any sequence matching a known
#'   allele (for either the given table or the current dataset) with a text
#'   name.
name_known_sequences <- function(results, known_alleles, name_args) {
  # Name all of the called alleles across samples
  results$summary <- name_alleles_in_table(results$summary, known_alleles,
                                           name_args)

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
#' @param results results list as produced by
#'   \code{\link{tidy_analyzed_dataset}}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{\link{load_locus_attrs}}.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.  \code{summary$Locus} is
#'   coerced to a factor with levels ordered according to their appearance in
#'   \code{locus_attrs$Locus}.  Order of rows in \code{summary} and entries in
#'   \code{data} are updated accordingly.
sort_results <- function(results, locus_attrs) {
  results$summary$Locus <- factor(results$summary$Locus,
                                  levels = locus_attrs$Locus)
  results$summary$Locus <- droplevels(results$summary$Locus)
  ord <- order_entries(results$summary)
  results$summary <- results$summary[ord, ]
  results$samples <- results$samples[ord]
  results
}
