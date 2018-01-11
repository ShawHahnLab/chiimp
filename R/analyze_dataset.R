# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze_sample and
# summarize_sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a list of processed samples
#' and a summary data frame showing each sample's summary row-by-row.
#'
#' @param dataset data frame of sample details as produced by
#'   \code{\link{prepare_dataset}}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{\link{load_locus_attrs}}.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match (see \code{\link{analyze_sample}}).
#' @param num.cores integer number of CPU cores to use in parallel for sample
#'   analysis.  Defaults to one less than half the number of detected cores with
#'   a minimum of 1.  If 1, the function will run without using the
#'   \code{parallel} package.
#' @param summary_args list of arguments to \code{summary.function}.
#' @param summary.function function to use when summarizing each sample's full
#'   details into the standard attributes Defaults to
#'   \code{\link{summarize_sample}}.
#'
#' @return list of results, with \code{summary} set to the single summary data
#'   frame and \code{data} the per-sample data frames.
#'
#' @export
analyze_dataset <- function(dataset,
                            locus_attrs,
                            nrepeats,
                            num.cores=max(1,
                                as.integer(parallel::detectCores() / 2) - 1),
                            summary_args,
                            summary.function=summarize_sample) {
  analyze.entry <- function(entry, summary.function) {
    seqs <- load_seqs(entry["Filename"])
    sample.data <- analyze_sample(seqs, locus_attrs, nrepeats)
    args <- c(list(sample.data = sample.data, locus.name = entry["Locus"]),
              summary_args)
    # If the dataset attributes specified expected lengths per sample per locus,
    # pass those in as well.
    # TODO just make the entry object a standard argument to these variant functions.
    expected_lengths <- unlist(unique(entry[c("ExpectedLength1",
                                              "ExpectedLength2")]))
    if (!is.null(expected_lengths) && !is.na(expected_lengths))
      args <- c(args, list(expected_lengths = expected_lengths))
    sample.summary <- do.call(summary.function, args)
    return(list(summary = sample.summary, data = sample.data))
  }
  if (num.cores > 1) {
    # Set up the cluster and export requried names.
    cluster_names <- c("locus_attrs",
                       "load_seqs",
                       "analyze_sample",
                       "find_matching_primer",
                       "check_motif",
                       "find_stutter",
                       "check_length")
    cluster <- parallel::makeCluster(num.cores)
    # https://stackoverflow.com/a/12232695/6073858
    parallel::clusterEvalQ(cluster, library(dnar))
    parallel::clusterExport(cluster, cluster_names)
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
  tidy_analyzed_dataset(dataset, raw.results)
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
  summaries <- lapply(raw.results, `[[`, 1)
  data <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind.data.frame, summaries)
  rownames(summary) <- rownames(dataset)
  names(data)       <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, data = data))
}
