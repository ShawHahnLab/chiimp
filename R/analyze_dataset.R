# Analyze all samples in a dataset.

# parallel execution of analyze.sample and summarize.sample
# binding sample summaries together in data frame

analyze.dataset <- function(dataset,
                            locus_attrs,
                            num.cores = max(1,
                                as.integer(parallel::detectCores() / 2) - 1)) {
  # Set up the cluster and export requried names.
  cluster <- parallel::makeCluster(num.cores)
  cluster_names <- c("locus_attrs",
         "load.seqs",
         "read.fastq",
         "read.fa",
         "summarize.sample",
         "analyze.sample",
         "find.matching.primer",
         "check.motif",
         "find.stutter",
         "check.length")
  parallel::clusterExport(cluster, cluster_names)
  # Load, analyze, and summarize each sample across the cluster.
  raw.results <- parallel::parApply(cluster, dataset, 1, function(entry) {
    seqs <- load.seqs(entry["Filename"])
    sample.data <- analyze.sample(seqs, locus_attrs, 3)
    sample.summary <- summarize.sample(sample.data, entry["Locus"])
    return(list(summary = sample.summary, data = sample.data))
  })
  parallel::stopCluster(cluster)
  # Recombined results into a summary data frame and a list of full sample data.
  tidy.analyzed.dataset(dataset, raw.results)
}

tidy.analyzed.dataset <- function(dataset, raw.results) {
  summaries <- lapply(raw.results, `[[`, 1)
  data <- lapply(raw.results, `[[`, 2)
  summary <- do.call(rbind.data.frame, summaries)
  rownames(summary) <- rownames(dataset)
  names(data)       <- rownames(dataset)
  return(list(summary = summary, data = data))
}
