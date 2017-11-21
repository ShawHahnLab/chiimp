# Analyze all samples in a dataset.

# These functions handle parallel execution of analyze_sample and
# summarize_sample and the binding sample summaries together in a data frame.

#' Analyze all samples in a dataset
#'
#' Load all samples for a given dataset and produce a list of processed samples
#' and a summary data frame showing each sample's summary row-by-row.
#'
#' @param dataset data frame of sample details as produced by
#'   \code{prepare_dataset}.
#' @param locus_attrs data frame of locus attributes as produced by
#'   \code{load_locus_attrs}.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match (see \code{analyze_sample}).
#' @param fraction.min numeric threshold for the minimum proportion of counts a
#'   given entry must have, compared to the total matching all criteria for that
#'   locus, to be considered as a potential allele.
#' @param counts.min numeric threshold for the minimum number of counts that
#'   must be present, in total across entries passing all filters, for potential
#'   alleles to be considered.
#' @param num.cores integer number of CPU cores to use in parallel for sample
#'   analysis.
#' @param summary.function function to use when summarizing each sample's full
#'   details into the standard attributes.
#'
#' @return list of results
#'
#' @export
analyze_dataset <- function(dataset,
                            locus_attrs,
                            nrepeats,
                            fraction.min,
                            counts.min,
                            num.cores=max(1,
                                as.integer(parallel::detectCores() / 2) - 1),
                            summary.function=summarize_sample) {
  analyze.entry <- function(entry, summary.function) {
    seqs <- load_seqs(entry["Filename"])
    sample.data <- analyze_sample(seqs, locus_attrs, nrepeats)
    sample.summary <- summary.function(sample.data,
                                       entry["Locus"],
                                       fraction.min,
                                       counts.min)
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
      # Load, analyze, and summarize each sample across the cluster.
      raw.results <- parallel::parApply(cluster, dataset, 1, analyze.entry,
                                        summary.function = summary.function)
    }, finally = {
      parallel::stopCluster(cluster)
    })
  } else {
    raw.results <- apply(dataset, 1, analyze.entry,
                         summary.function = summary.function)
  }
  # Recombined results into a summary data frame and a list of full sample data.
  tidy_analyzed_dataset(dataset, raw.results)
}

# rearrange the pairs of sample summary / sample data objects into a single
# summary cross-sample data frame and a list of detailed per-sample data frames.
tidy_analyzed_dataset <- function(dataset, raw.results) {
  summaries <- lapply(raw.results, `[[`, 1)
  data <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind.data.frame, summaries)
  rownames(summary) <- rownames(dataset)
  names(data)       <- rownames(dataset)
  summary <- cbind(dataset, summary)
  return(list(summary = summary, data = data))
}
