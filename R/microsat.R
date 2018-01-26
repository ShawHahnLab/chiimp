# Run a full microsatellite analysis and handle configuration and command-line
# execution.

#' Analyze Microsatellites
#'
#' Analyze DNA microsatellites in high-throughput sequencing datasets.
"_PACKAGE"


# Configuration -----------------------------------------------------------

#' Default microsatellite configuration
#'
#' The entries in this list show the build-time configuration defaults for all
#' aspects of the microsatellite analysis.  These can be overridden by passing a
#' list to \code{\link{full_analysis}} with entries of the same names.
#'
#' Notable Options:
#'   * dataset_opts:
#'     * dp: directory path to input sequence files
#'     * pattern: regular expression for the input filename pattern
#'     * ord: order of fields in the input filename pattern
#'   * output:
#'     * dp: directory path for saving output data
#'   * fp_locus_attrs: file path to locus attributes TSV file
#'   * fp_genotypes_known: file path to known genotypes TSV file
#' @md
#'
#' @export
config.defaults <- list(
  ## Arguments to prepare_dataset
  dataset_opts = list(
    dp = "str-data",
    pattern = "(\\d+)-(\\d+)-([A-Za-z0-9]+).fast[aq](?:\\.gz)",
    ord = c(1, 2, 3),
    autorep = FALSE
  ),
  # Other input and output paths
  fp_locus_attrs = "locus_attrs.tsv",
  fp_allele_names = NULL,
  fp_genotypes_known = NULL,
  ## Names for output files and directories
  output = list(
    dp = "str-results",  # Main output directory
    fp_summary = "summary.csv",  # Dataset summary table
    fp_report = "report.html",  # Report document
    fp_dist_mat = "sample-distances.csv",  # Sample-to-sample distances
    fp_rds = NULL,  # Data file to save results to
    dp_histograms = "histograms",  # Read count by length histograms
    dp_alignments = "alignments",  # Sequence alignments across alleles
    dp_alignment_images = "alignment-images",  # Images of alignments
    dp_processed_samples = "processed-samples",  # Sample data tables
    dp_allele_seqs = "allele-sequences"  # FASTA sequences for alleles
  ),
  ## Options for analyze_dataset
  dataset_analysis = list(ncores = 0),
  ## Sample genotyping settings
  sample_analysis = list(nrepeats = 3),
  sample_summary_func = "summarize_sample",
  sample_summary = list(fraction.min = 0.05,
                        counts.min = 500),
  ## Report generation settings
  # Should a report be generated?
  report = TRUE,
  # Should the code executed in the report generation be included in the report?
  report.echo = FALSE,
  # Title and other metadata at top of the report.
  report.title = "Microsatellite Report",
  report.author = NULL,
  # Length of suffix on automated allele names in tables.
  report.hash_len = 6,
  # List of vectors of locus names to use to break up tables into reasonable
  # sizes and/or order locus names explicitly.
  report.locus_chunks = NULL,
  # Group together genotype table rows by sample?
  report.group_samples = FALSE,
  # Text to use for NA entries in Replicates column of tables (i.e. Pooled)
  report.na.replicates = "",
  # Parameters controlling how identifications are reported: dist_range is how
  # closeby (to the closest case) the next-nearest individuals must be to be
  # listed as similar to a sample, and dist_max is the maximum distance for a
  # given individual to be listed.
  report.dist_range = 2,
  report.dist_max = 3,
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
  verbose = TRUE)

# Analysis ----------------------------------------------------------------

#' Perform a full microsatellite analysis
#'
#' Given a list of configuration options, run all aspects of a microsatellite
#' analysis, and save the corresponding output files.
#'
#' @param config list of configuration options.  See the summary for
#'   \code{\link{config.defaults}} for more details.
#' @param dataset optional custom data frame to override the sample attributes
#'   defined in the config.
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @export
full_analysis <- function(config, dataset=NULL) {
  # Overaly explicit configuration onto the default settings
  config_full <- utils::modifyList(config.defaults, config)

  # Make output path absolute
  config_full$output$dp <-
    if (substr(config_full$output$dp, 1, 1) != .Platform$file.sep) {
      file.path(normalizePath("."), config_full$output$dp)
    } else {
      config_full$output$dp
    }

  # Only show identifications if a known genotypes table was supplied
  config_full$report.sections$identifications <-
    ! is.null(config_full$fp_genotypes_known) &&
    config_full$report.sections$identifications

  cfg <- config_full
  if (cfg$verbose)
    logmsg(paste0("Loading dataset: ", cfg$dataset_opts$dp, "..."))
  if (is.null(dataset)) {
    dataset <- do.call(prepare_dataset, cfg$dataset_opts)
  }
  if (cfg$verbose)
    logmsg(paste0("Loading locus attrs: ", cfg$fp_locus_attrs, "..."))
  locus_attrs <- load_locus_attrs(cfg$fp_locus_attrs)
  if (cfg$verbose) logmsg("Analyzing samples...")
  idx <- match(cfg$sample_summary_func, sample_summary_funcs, nomatch = 1)
  sample_summary_func <- get(sample_summary_funcs[idx])
  results <- analyze_dataset(dataset, locus_attrs,
                             nrepeats = cfg$sample_analysis$nrepeats,
                             ncores = cfg$dataset_analysis$ncores,
                             summary_args = cfg$sample_summary,
                             summary.function = sample_summary_func)
  # Reorder entries and levels to match locus_attrs.
  # TODO merge these steps into analyze_dataset or summarize_dataset
  results$summary$Locus <- factor(results$summary$Locus,
                                  levels = rownames(locus_attrs))
  ord <- order_entries(results$summary)
  results$summary <- results$summary[ord, ]
  results$data <- results$data[ord]
  results$locus_attrs <- locus_attrs
  if (cfg$verbose) logmsg("Summarizing results...")
  genotypes.known <- NULL
  if (!is.null(cfg$fp_genotypes_known))
    genotypes.known <- load_genotypes(cfg$fp_genotypes_known)
  results <- summarize_dataset(results, genotypes.known)
  if (!is.null(cfg$fp_genotypes_known))
    results$closest_matches <- find_closest_matches(results$dist_mat_known,
                                                range = cfg$report.dist_range,
                                                maximum = cfg$report.dist_max)
  results$cts_per_locus <- tally_cts_per_locus(results)
  results$config <- config_full
  results$allele.names <- NULL
  if (!is.null(cfg$fp_allele_names))
    results$allele.names <- load_allele_names(cfg$fp_allele_names)
  if (cfg$verbose) logmsg("Saving output files...")
    save_data(results, results$config)
  if (cfg$report) {
    if (cfg$verbose) logmsg("Creating report...")
    render_report(results, results$config)
  }
  if (cfg$verbose) logmsg("Done.")
  return(results)
}

#' Handle full microsatellite analysis from command-line
#'
#' Read configuration from command-line arguments and run
#' \code{\link{full_analysis}}.
#'
#' @param args optional character vector of arguments to use rather than those
#'   detected with \code{\link[base]{commandArgs}}.
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @export
main <- function(args=NULL) {
  if (missing(args))
    args <- commandArgs(trailingOnly = TRUE)
  desc <- "Identify microsatellite alleles fasta/fastq files"
  p <- argparser::arg_parser(desc)
  p <- argparser::add_argument(p, "config", help = "configuration file path")
  args_parsed <- argparser::parse_args(p, args)
  config <- load_config(args_parsed$config)
  full_analysis(config)
}

# Util --------------------------------------------------------------------

#' Render Microsatellite Report Document
#'
#' Write an HTML report summarizing the results of a full microsatellite
#' analysis.
#'
#' @param results list of microsatellite analysis results as produced by
#'   \code{\link{full_analysis}}.
#' @param config list of configuration options (see
#'   \code{\link{config.defaults}}).
render_report <- function(results, config) {
  with(config, {
    fp_report_in <- system.file("report", "report.Rmd", package = "microsat")
    fp_report_out <- file.path(output$dp, output$fp_report)
    if (!dir.exists(dirname(fp_report_out)))
      dir.create(dirname(fp_report_out), recursive = TRUE)
    pandoc_metadata <- c(title = report.title,
                         author = report.author,
                         date = format(Sys.Date(), "%Y-%m-%d"))
    pandoc_args <- format_pandoc_args(pandoc_metadata)
    pandoc_args <- c(pandoc_args, paste0("--css=", "report.css"))
    rmarkdown::render(fp_report_in, quiet = TRUE, output_file = fp_report_out,
                      output_options = list(pandoc_args = pandoc_args))
  })
}

#' Save Microsatellite Analysis to Disk
#'
#' Save the various parts of a finished microsatellite analysis to data files
#' and plot images.
#'
#' @param results list of microsatellite analysis results as produced by
#'   \code{\link{full_analysis}}.
#' @param config list of configuration options (see
#'   \code{\link{config.defaults}}).
save_data <- function(results, config) {
  save_histograms(results,
              file.path(config$output$dp, config$output$dp_histograms))
  save_results_summary(results$summary,
              file.path(config$output$dp, config$output$fp_summary))
  save_alignments(results$alignments,
              file.path(config$output$dp, config$output$dp_alignments))
  save_alignment_images(results$alignments,
              file.path(config$output$dp, config$output$dp_alignment_images))
  save_sample_data(results$data,
              file.path(config$output$dp, config$output$dp_processed_samples))
  save_allele_seqs(results$summary,
              file.path(config$output$dp, config$output$dp_allele_seqs))
  save_dist_mat(results$dist_mat,
              file.path(config$output$dp, config$output$fp_dist_mat))
  if (!is.null(config$output$fp_rds))
    saveRDS(results,
              file.path(config$output$dp, config$output$fp_rds))
}

#' Create Pandoc Metadata Argument Strings
#'
#' Convert named list of pandoc metadata options into command-line
#' \code{--metadata} argument strings.
#'
#' @param metadata named list of
#'   \href{http://pandoc.org/MANUAL.html#metadata-blocks}{pandoc metadata
#'   options}.
#'
#' @return list of \code{--metadata} argument strings
format_pandoc_args <- function(metadata) {
  metadata <- paste(names(metadata),
                    lapply(metadata,
                           function(s) paste0("\"", s, "\"")), sep = ":")
  paste("--metadata=", metadata, sep = "")
}

#' Write Log Message
#'
#' Print a log message to the standard error stream.
#'
#' @param msg text to print.
#' @param col2 extra text to show at right margin; defaults to current time.
#' @param end ending to concatenate to message; defaults to newline character.
logmsg <- function(msg, col2=as.character(Sys.time()), end="\n") {
  if (!is.null(col2)) {
    # right-justify col2, aim to fit total msg within 80 characters
    pad <- max(1, 80 - nchar(msg) - nchar(col2))
    msg <- paste0(msg, paste(rep(" ", pad),  collapse = ""), col2)
  }
  # stderr: file descriptor 2
  cat(paste0(msg, end), file = 2)
}
