# Run a full microsatellite analysis and handle configuration and command-line
# execution.

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
  dp.output="str-results",
  fp.locus_attrs="locus_attrs.tsv",
  fp.output.summary="summary.csv",
  dp.output.alignments="alignments",
  dp.output.alignment_images="alignment-images",
  dp.output.processed_samples="processed-samples")

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
  config_full <- modifyList(config.defaults, config)
  with(config_full, {
    dataset <- prepare.dataset(dp.data, pattern)
    locus_attrs <- load.locus_attrs(fp.locus_attrs)
    results <- analyze.dataset(dataset, locus_attrs)
    results <- summarize.dataset(results)
    save.results_summary(results$summary, file.path(dp.output, fp.output.summary))
    save.alignments(results$alignments, file.path(dp.output, dp.output.alignments))
    save.alignment_images(results$alignments, file.path(dp.output, dp.output.alignment_images))
    save.all_sample_data(results$data, file.path(dp.output, dp.output.processed_samples))
    results$config <- config_full
    # TODO histograms
    # TODO locus performance
    # TODO alleles as fasta files
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
  p <- argparser::add_argument("config", help = "configuration file path")
  args_parsed <- argparser::parse_args(p, args)
  config <- load.config(args_parsed$config)
  full_analysis(config)
}
