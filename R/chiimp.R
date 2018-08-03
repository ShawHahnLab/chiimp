# Run a full microsatellite analysis and handle configuration and command-line
# execution.

#' Analyze Microsatellites
#'
#' Analyze DNA microsatellites in high-throughput sequencing datasets.
#'
#' @details
#'
#' Starting from file inputs and producing file outputs, the overall workflow
#' (handled by \code{\link{full_analysis}} as a configuration-driven wrapper for
#' the entire process) is:
#'
#' * Load input data
#'   * Load data frame of sample information from a spreadsheet via
#'     \code{\link{load_dataset}} or directly from filenames via
#'     \code{\link{prepare_dataset}}.
#'   * Load data frame of locus attributes via \code{\link{load_locus_attrs}}
#'   * Optionally, load data frame of names for allele sequences via
#'     \code{\link{load_allele_names}}.
#'   * Optionally, load data frame of known genotypes for named individuals via
#'     \code{\link{load_genotypes}}.
#' * Analyze dataset via \code{\link{analyze_dataset}}
#'   * Load each sequence data file into a character vector with
#'     \code{\link{load_seqs}} and process into a dereplicated data frame with
#'     \code{\link{analyze_seqs}}.
#'   * For each sample, filter the sequences from the relevant per-file data
#'     frame to just those matching the expected locus and identify possible
#'     alleles, via \code{\link{analyze_sample}}.  (There may be a many-to-one
#'     relationship of samples to files, for example with sequencer
#'     multiplexing.)
#'   * Process each per-sample data frames into a summary list of attributes
#'     giving alleles identified and related information, via
#'     \code{\link{summarize_sample}}.
#'   * Organize \code{analyze_dataset} results into a list of per-file data
#'     frames, a list of per-sample data frames, and a single summary data
#'     frame across all samples.
#' * Summarize results and add additional comparisons (cross-sample and to
#'   known-individual) via \code{\link{summarize_dataset}}.
#'   * Tabulate sequence counts per sample matching each locus' primer via
#'     \code{\link{tally_cts_per_locus}}.
#'   * Align identified alleles for each locus via
#'     \code{\link{align_alleles}}.
#'   * Create a sample-to-sample distance matrix of allele mismatches via
#'     \code{\link{make_dist_mat}}.
#'   * If genotypes for known individuals were provided, create a
#'     sample-to-known-individual distance matrix via
#'     \code{\link{make_dist_mat_known}}.
#'   * If identities of samples were provided, score genotyping success via
#'     \code{\link{match_known_genotypes}} and
#'     \code{\link{categorize_genotype_results}}.
#' * Save analysis results to files.
#' * Create an HTML report document summarizing all results.
#'
#' The workflow above outlines CHIIMP's behavior when called as a standalone
#' program, where \code{\link{main}} loads a configuration file into a nested
#' list of options and calls \code{\link{full_analysis}}.  The public functions
#' linked above can also be used idependently; see the documentation and code
#' examples for the individual functions for more information.
#'
#' @md
#'
"_PACKAGE"

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
    logmsg(paste0("Loading dataset: ", ifelse(is.null(cfg$fp_dataset),
                                              cfg$dataset_opts$dp,
                                              cfg$fp_dataset), "..."))
  if (is.null(dataset)) {
    if (is.null(cfg$fp_dataset)) {
      dataset <- do.call(prepare_dataset, cfg$dataset_opts)
    } else {
      dataset <- load_dataset(cfg$fp_dataset)
    }
  }

  if (cfg$verbose)
    logmsg(paste0("Loading locus attrs: ", cfg$fp_locus_attrs, "..."))
  locus_attrs <- load_locus_attrs(cfg$fp_locus_attrs)

  if (cfg$verbose) logmsg("Analyzing samples...")
  idx <- match(cfg$sample_summary_func, sample_summary_funcs, nomatch = 1)
  sample_summary_func <- get(sample_summary_funcs[idx])
  idx <- match(cfg$sample_analysis_func, sample_analysis_funcs, nomatch = 1)
  sample_analysis_func <- get(sample_analysis_funcs[idx])
  allele.names <- NULL
  if (!is.null(cfg$fp_allele_names))
    allele.names <- load_allele_names(cfg$fp_allele_names)
  results <- analyze_dataset(dataset, locus_attrs,
                             nrepeats = cfg$seq_analysis$nrepeats,
                             ncores = cfg$dataset_analysis$ncores,
                             analysis_opts = cfg$sample_analysis_opts,
                             summary_opts = cfg$sample_summary_opts,
                             analysis_function = sample_analysis_func,
                             summary_function = sample_summary_func,
                             known_alleles = allele.names,
                             name_args = cfg$dataset_analysis$name_args)
  results$allele.names <- allele.names
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
  results$config <- config_full
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
    fp_report_in <- system.file("report", "report.Rmd", package = "chiimp")
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
  save_seqfile_data(results$files,
              file.path(config$output$dp, config$output$dp_processed_files))
  save_sample_data(results$samples,
                    file.path(config$output$dp,
                              config$output$dp_processed_samples))
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
