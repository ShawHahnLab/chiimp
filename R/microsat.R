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
#' @md
#'
#' @export
config.defaults <- list(
  dp.data="str-data",
  pattern="(\\d+)-(\\d+)-([A-Za-z0-9]+).fast[aq](?:\\.gz)",
  ord=c(1, 2, 3),
  fp.locus_attrs="locus_attrs.tsv",
  dp.output="str-results",
  # Names for files and subdirectories under dp.output
  fp.output.summary="summary.csv",
  fp.report="report.html",
  fp.output.dist_mat="sample-distances.csv",
  dp.output.alignments="alignments",
  dp.output.alignment_images="alignment-images",
  dp.output.processed_samples="processed-samples",
  dp.output.allele_seqs="allele-sequences",
  # Report generation settings
  report=TRUE,
  report.echo=FALSE,
  report.sections = c(genotypes  = TRUE,
                      distances  = TRUE,
                      flags      = TRUE,
                      alignments = TRUE),
  # Other settings
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
  config.full <- modifyList(config.defaults, config)

  # Make output path absolute
  config.full$dp.output <-
    if (substr(config.full$dp.output, 1, 1) != .Platform$file.sep) {
      file.path(normalizePath('.'), config.full$dp.output)
    } else {
      config.full$dp.output
    }

  with(config.full, {
    if (verbose) logmsg(paste("Loading dataset from", dp.data, "..."))
    dataset <- prepare_dataset(dp.data, pattern)
    if (verbose)
      logmsg(paste("Loading locus attributes from", fp.locus_attrs, "..."))
    locus_attrs <- load_locus_attrs(fp.locus_attrs)
    if (verbose) logmsg("Analyzing samples...")
    results <- analyze_dataset(dataset, locus_attrs)
    if (verbose) logmsg("Summarizing results...")
    results <- summarize_dataset(results)
    results$config <- config.full
    if (verbose) logmsg("Saving output files...")
    save_results_summary(results$summary, file.path(dp.output, fp.output.summary))
    save_alignments(results$alignments, file.path(dp.output, dp.output.alignments))
    save_alignment_images(results$alignments, file.path(dp.output, dp.output.alignment_images))
    save_sample_data(results$data, file.path(dp.output, dp.output.processed_samples))
    save_allele_seqs(results$summary, file.path(dp.output, dp.output.allele_seqs))
    save_dist_mat(results$dist_mat, file.path(dp.output, fp.output.dist_mat))
    if (report) {
      if (verbose) logmsg("Creating report...")
      fp.report.in <- system.file("report", "report.Rmd", package="microsat")
      fp.report.out <- file.path(dp.output, fp.report)
      if (!dir.exists(dirname(fp.report.out)))
        dir.create(dirname(fp.report.out), recursive = TRUE)
      rmarkdown::render(fp.report.in, quiet = TRUE, output_file = fp.report.out)
    }
    # TODO histograms
    # TODO locus performance
    if (verbose) logmsg("Done.")
    return(results)
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
