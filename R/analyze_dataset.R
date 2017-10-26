# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze.sample and
# summarize.sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a list of processed samples
#' and a summary data frame showing each sample's summary row-by-row.
#'
#' @param dataset data frame of sample details as produced by
#'   \code{prepare.dataset}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{load.locus_attrs}.
#' @param num.cores integer number of CPU cores to use in parallel for sample
#'   analysis.
#' @param summary.function function to use when summarizing each sample's full
#'   details into the standard attributes.
#'
#' @return list of results
#'
#' @export
analyze.dataset <- function(dataset,
                            locus_attrs,
                            num.cores=max(1,
                                as.integer(parallel::detectCores() / 2) - 1),
                            summary.function=summarize.sample) {
  analyze.entry <- function(entry) {
    seqs <- load.seqs(entry["Filename"])
    sample.data <- analyze.sample(seqs, locus_attrs, 3)
    sample.summary <- summary.function(sample.data, entry["Locus"])
    return(list(summary = sample.summary, data = sample.data))
  }
  if (num.cores > 1) {
    # Set up the cluster and export requried names.
    cluster <- parallel::makeCluster(num.cores)
    cluster_names <- c("locus_attrs",
                       "load.seqs",
                       "summarize.sample",
                       "summarize.sample.by_length",
                       "analyze.sample",
                       "find.matching.primer",
                       "check.motif",
                       "find.stutter",
                       "check.length")
    # https://stackoverflow.com/a/12232695/6073858
    parallel::clusterEvalQ(cluster, library(dnar))
    parallel::clusterExport(cluster, cluster_names)
    # Load, analyze, and summarize each sample across the cluster.
    raw.results <- parallel::parApply(cluster, dataset, 1, analyze.entry)
    parallel::stopCluster(cluster)
  } else {
    raw.results <- apply(dataset, 1, analyze.entry)
  }
  # Recombined results into a summary data frame and a list of full sample data.
  tidy.analyzed.dataset(dataset, raw.results)
}

tidy.analyzed.dataset <- function(dataset, raw.results) {
  summaries <- lapply(raw.results, `[[`, 1)
  data <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind.data.frame, summaries)
  rownames(summary) <- rownames(dataset)
  names(data)       <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, data = data))
}
