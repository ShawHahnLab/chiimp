# Run a full microsatellite analysis and handle configuration and command-line
# execution.

#' Analyze Microsatellites
#'
#' Analyze DNA microsatellites in high-throughput sequencing datasets.
"_PACKAGE"

#' Default microsatellite configuration
#'
#' The entries in this list show the build-time configuration defaults for all
#' aspects of the microsatellite analysis.  These can be overridden by passing a
#' list to \code{full_analysis} with entries of the same names.
#'
#' Notable Options:
#'   * dp.data: directory path to input sequence files
#'   * pattern: regular expression for the input filename pattern
#'   * ord: order of fields in the input filename pattern
#'   * dp.output: directory path for saving output data
#'   * fp.locus_attrs: file path to locus attributes TSV file
#'   * fp.genotypes.known: file path to known genotypes TSV file
#' @md
#'
#' @export
config.defaults <- list(
  dp.data="str-data",
  pattern="(\\d+)-(\\d+)-([A-Za-z0-9]+).fast[aq](?:\\.gz)",
  ord=c(1, 2, 3),
  fp.locus_attrs="locus_attrs.tsv",
  fp.allele.names=NULL,
  fp.genotypes.known=NULL,
  dp.output="str-results",
  ## Names for files and subdirectories under dp.output
  fp.output.summary="summary.csv",  # Dataset summary table
  fp.report="report.html",  # Report document
  fp.output.dist_mat="sample-distances.csv",  # Sample-to-sample distances
  fp.output.rds=NULL,  # Data file to save results to
  dp.output.histograms="histograms",  # Read count by length histograms
  dp.output.alignments="alignments",  # Sequence alignments across alleles
  dp.output.alignment_images="alignment-images",  # Visualizations of alignments
  dp.output.processed_samples="processed-samples",  # Sample data tables
  dp.output.allele_seqs="allele-sequences",  # FASTA sequences for alleles
  ## Sample genotyping settings
  sample_analysis = list(nrepeats = 3),
  sample_summary_func = "summarize_sample",
  sample_summary = list(fraction.min = 0.05,
                        counts.min = 500),
  ## Report generation settings
  # Should a report be generated?
  report=TRUE,
  # Should the code executed in the report generation be included in the report?
  report.echo=FALSE,
  # Title at top of the report.
  report.title="Microsatellite Report",
  # Length of suffix on automated allele names in tables.
  report.hash_len=6,
  # List of vectors of locus names to use to break up tables into reasonable
  # sizes and/or order locus names explicitly.
  report.locus_chunks=NULL,
  # Group together genotype table rows by sample?
  report.group_samples=FALSE,
  # Text to use for NA entries in Replicates column of tables (i.e. Pooled)
  report.na.replicates="",
  # Parameters controlling how identifications are reported: dist_range is how
  # closeby (to the closest case) the next-nearest individuals must be to be
  # listed as similar to a sample, and dist_max is the maximum distance for a
  # given individual to be listed.
  report.dist_range=2,
  report.dist_max=8,
  # Sub-sections of the report that can be excluded by overriding these with
  # FALSE.
  report.sections = list(genotypes       = TRUE,
                         identifications = TRUE,
                         distances       = TRUE,
                         flags           = TRUE,
                         alignments      = TRUE,
                         contamination   = TRUE),
  ## Other settings
  # Print status messages to standard error.
  verbose=TRUE)

#' Perform a full microsatellite analysis
#'
#' Given a list of configuration options, run all aspects of a microsatellite
#' analysis, and save the corresponding output files.
#'
#' @param config list of configuration options.  See the summary for
#'   \code{config.defaults} for more details.
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @export
full_analysis <- function(config) {
  # Overaly explicit configuration onto the default settings
  config.full <- utils::modifyList(config.defaults, config)

  # Make output path absolute
  config.full$dp.output <-
    if (substr(config.full$dp.output, 1, 1) != .Platform$file.sep) {
      file.path(normalizePath('.'), config.full$dp.output)
    } else {
      config.full$dp.output
    }

  # Only show identifications if a known genotypes table was supplied
  config.full$report.sections$identifications <-
    ! is.null(config.full$fp.genotypes.known) &&
    config.full$report.sections$identifications

  with(config.full, {
    if (verbose) logmsg(paste0("Loading dataset from ", dp.data, "..."))
    dataset <- prepare_dataset(dp.data, pattern, ord)
    if (verbose)
      logmsg(paste0("Loading locus attributes from ", fp.locus_attrs, "..."))
    locus_attrs <- load_locus_attrs(fp.locus_attrs)
    if (verbose) logmsg("Analyzing samples...")
    idx <- match(sample_summary_func, sample_summary_funcs, nomatch = 1)
    sample_summary_func <- get(sample_summary_funcs[idx])
    results <- analyze_dataset(dataset, locus_attrs,
                               nrepeats = sample_analysis$nrepeats,
                               fraction.min = sample_summary$fraction.min,
                               counts.min = sample_summary$counts.min,
                               summary.function = sample_summary_func)
    # Reorder entries and levels to match locus_attrs.
    ord <- order(match(dataset$Locus, rownames(locus_attrs)),
                 order_entries(dataset))
    results$summary <- results$summary[ord, ]
    results$data <- results$data[ord]
    results$locus_attrs <- locus_attrs
    levels(results$summary$Locus) <- rownames(locus_attrs)
    if (verbose) logmsg("Summarizing results...")
    genotypes.known <- NULL
    if (!is.null(fp.genotypes.known))
      genotypes.known <- load_genotypes(fp.genotypes.known)
    results <- summarize_dataset(results, genotypes.known)
    results$config <- config.full
    if (verbose) logmsg("Saving output files...")
    save_histograms(results, file.path(dp.output, dp.output.histograms))
    save_results_summary(results$summary, file.path(dp.output, fp.output.summary))
    save_alignments(results$alignments, file.path(dp.output, dp.output.alignments))
    save_alignment_images(results$alignments, file.path(dp.output, dp.output.alignment_images))
    save_sample_data(results$data, file.path(dp.output, dp.output.processed_samples))
    save_allele_seqs(results$summary, file.path(dp.output, dp.output.allele_seqs))
    save_dist_mat(results$dist_mat, file.path(dp.output, fp.output.dist_mat))
    if (!is.null(fp.output.rds))
      saveRDS(results, file.path(dp.output, fp.output.rds))
    if (report) {
      if (verbose) logmsg("Creating report...")
      render_report(results, results$config)
    }
    if (verbose) logmsg("Done.")
    return(results)
  })
}

render_report <- function(results, config) {
  with(config, {
    fp.report.in <- system.file("report", "report.Rmd", package="microsat")
    fp.report.out <- file.path(dp.output, fp.report)
    if (!dir.exists(dirname(fp.report.out)))
      dir.create(dirname(fp.report.out), recursive = TRUE)
    pandoc_args <- c(paste0("--metadata=title:\"", report.title, "\""))
    rmarkdown::render(fp.report.in, quiet = TRUE, output_file = fp.report.out,
                      output_options = list(pandoc_args = pandoc_args))
  })
}

#' Handle full microsatellite analysis from command-line
#'
#' Read configuration from command-line arguments and run \code{full_analysis}.
#'
#' @param args optional character vector of arguments to use rather than those
#'   detected with \code{commandArgs}.
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @export
main <- function(args=NULL) {
  if (missing(args))
    args <- commandArgs(trailingOnly = TRUE)
  p <- argparser::arg_parser("Identify microsatellite alleles fasta/fastq files")
  p <- argparser::add_argument(p, "config", help = "configuration file path")
  args_parsed <- argparser::parse_args(p, args)
  config <- load_config(args_parsed$config)
  full_analysis(config)
}

logmsg <- function(msg) {
  cat(paste0(msg, "\n"), file="/dev/stderr")
}
