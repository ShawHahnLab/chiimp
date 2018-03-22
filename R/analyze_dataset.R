# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze_sample and
# summarize_sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a list of processed samples
#' and a summary data frame showing each sample's summary row-by-row.  The
#' entries in the processed-samples list and the rows in the summary data frame
#' will be sorted according to the ordering of loci in \code{locus_attrs} and
#' by the sample attributes.
#'
#' @param dataset data frame of sample details as produced by
#'   \code{\link{prepare_dataset}}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{\link{load_locus_attrs}}.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match (see \code{\link{analyze_sample}}).
#' @param ncores integer number of CPU cores to use in parallel for sample
#'   analysis.  Defaults to one less than half the number of detected cores with
#'   a minimum of 1.  If 1, the function will run without using the
#'   \code{parallel} package.
#' @param summary_args list of supplemental arguments to
#'   \code{summary.function}.
#' @param summary.function function to use when summarizing each sample's full
#'   details into the standard attributes.  Defaults to
#'   \code{\link{summarize_sample}}.
#' @param known_alleles data frame of custom allele names as defined for
#'   \code{\link{load_allele_names}}.  if NULL only the names automatically
#'   generated for the dataset summary will be used.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.
#'
#' @export
analyze_dataset <- function(dataset,
                            locus_attrs,
                            nrepeats,
                            ncores = 0,
                            summary_args,
                            summary.function=summarize_sample,
                            known_alleles=NULL) {
  if (ncores == 0) {
    ncores <- max(1, as.integer(parallel::detectCores() / 2) - 1)
  }
  analyze.entry <- function(entry, summary.function) {
    seqs <- load_seqs(entry["Filename"])
    sample.data <- analyze_sample(seqs, locus_attrs, nrepeats)
    args <- c(list(sample.data = sample.data, sample.attrs = entry),
              summary_args)
    sample.summary <- do.call(summary.function, args)
    return(list(summary = sample.summary, data = sample.data))
  }
  if (ncores > 1) {
    # Set up the cluster and export required names.
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
    parallel::clusterExport(cl = cluster,
                            varlist = cluster_names,
                            envir = environment())
    tryCatch({
      # Load, analyze, and summarize each sample across the cluster.  Each row
      # in the dataset data frame will be given as the entry argument to
      # analyze.entry.
      raw.results <- parallel::parApply(cluster, dataset, 1, analyze.entry,
                                        summary.function = summary.function)
    },
    finally = {
      parallel::stopCluster(cluster)
    })
  } else {
    raw.results <- apply(dataset, 1, analyze.entry,
                         summary.function = summary.function)
  }
  # Recombined results into a summary data frame and a list of full sample data.
  results <- tidy_analyzed_dataset(dataset, raw.results)
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
#' @param dataset data frame of sample details as produced by
#'   \code{\link{prepare_dataset}}.
#' @param raw.results list of pairs of sample summary and sample data data
#'   frames (from \code{\link{summarize_sample}} and
#'   \code{\link{analyze_sample}}).
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.
tidy_analyzed_dataset <- function(dataset, raw.results) {
  summaries <- lapply(raw.results, function(s) {
    data.frame(s[[1]], stringsAsFactors = FALSE)
  })
  data <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind, summaries)
  rownames(summary) <- rownames(dataset)
  names(data)       <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, data = data))
}

#' Name known allele sequences
#'
#' For the given results list (pair of summary data frame and list of per-sample
#' data frames as produced by \code{\link{tidy_analyzed_dataset}}), add columns
#' to all data frames defining names for recognized sequences.  For the summary
#' data frame this will be Allele1Name and Allele2Name.  For each sample data
#' frame this will be SeqName, defined for any sequences represented in the
#' summary or in a given known alleles set.
#'
#' @param results results list as produced by
#'   \code{\link{tidy_analyzed_dataset}}.
#' @param known_alleles data frame of custom allele names as defined for
#'   \code{\link{load_allele_names}}.  if NULL only the names automatically
#'   generated for the summary will be used.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.  A "SeqName" column in
#'   sample data frames and "Allele1Name" and "Allele2Name" columns in the
#'   summary data frame will associate any sequence matching a known allele (for
#'   either the given table or the current dataset) with a text name.
name_known_sequences <- function(results, known_alleles) {
  # Name all of the called alleles across samples
  results$summary <- name_alleles_in_table(results$summary, known_alleles)

  # Create table of allele names for current dataset
  a1 <- results$summary[, c("Locus", "Allele1Seq", "Allele1Name")]
  a2 <- results$summary[, c("Locus", "Allele2Seq", "Allele2Name")]
  colnames(a1) <- c("Locus", "Seq", "Name")
  colnames(a2) <- c("Locus", "Seq", "Name")
  aT <- rbind(a1, a2)
  aT <- unique(aT[! is.na(aT$Seq), ])

  # Merge into given known alleles table (if present)
  known_alleles <- if (is.null(known_alleles)) {
    aT
  } else {
    known_alleles$Name <- as.character(known_alleles$Name)
    unique(rbind(known_alleles, aT))
  }

  # Name recognized sequences in each sample data frame
  results$data <- lapply(results$data, function(d) {
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
  results$data <- results$data[ord]
  results
}
