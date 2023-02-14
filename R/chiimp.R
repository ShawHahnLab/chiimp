# Run a full microsatellite analysis and handle configuration and command-line
# execution.

#' CHIIMP
#'
#' Computational, High-throughput Individual Identification through
#' Microsatellite Profiling.  For a conceptual overview see the latest
#' [user guide](https://shawhahnlab.github.io/chiimp/GUIDE.pdf) and
#' [additional documentation](https://shawhahnlab.github.io/chiimp/docs/) at
#' <https://shawhahnlab.github.io/chiimp/>.
#'
#' @details
#'
#' Starting from file inputs and producing file outputs, the overall workflow
#' (handled by [full_analysis] as a configuration-driven wrapper for
#' the entire process) is:
#'
#' * Load input data.  The input spreadsheets are text files using
#'   comma-separated values (CSV).
#'   * Load data frame of sample information from a spreadsheet via
#'   [load_dataset] or directly from filenames via [prepare_dataset].
#'   * Load data frame of locus attributes via [load_locus_attrs]
#'   * Optionally, load data frame of names for allele sequences via
#'   [load_allele_names].
#'   * Optionally, load data frame of known genotypes for named individuals via
#'     [load_genotypes].
#' * Analyze dataset via [analyze_dataset]
#'   * Load each sequence data file into a character vector with
#'     [load_seqs] and process into a dereplicated data frame with
#'     [analyze_seqs].
#'   * For each sample, filter the sequences from the relevant per-file data
#'     frame to just those matching the expected locus and identify possible
#'     alleles, via [analyze_sample].  (There may be a many-to-one relationship
#'     of samples to files, for example with sequencer multiplexing.)
#'   * Process each per-sample data frames into a summary list of attributes
#'     giving alleles identified and related information, via
#'     [summarize_sample].
#'   * Organize [analyze_dataset] results into a list of per-file data
#'     frames, a list of per-sample data frames, and a single summary data
#'     frame across all samples.
#' * Summarize results and add additional comparisons (cross-sample and to
#'   known-individual) via [summarize_dataset].
#'   * Tabulate sequence counts per sample matching each locus' primer via
#'     [tally_cts_per_locus].
#'   * Align identified alleles for each locus via [align_alleles].
#'   * Create a sample-to-sample distance matrix of allele mismatches via
#'     [make_dist_mat].
#'   * If genotypes for known individuals were provided, create a
#'     sample-to-known-individual distance matrix via [make_dist_mat_known].
#'   * If identities of samples were provided, score genotyping success via
#'     [match_known_genotypes] and [categorize_genotype_results].
#' * Save analysis results to files.  Spreadsheets are in CSV format for output
#'   as well as input.  Some output files are in FASTA format (alignments and
#'   alleles) or are PNG images (alignment visualization and sequence count
#'   histograms).  If specified in the configuration, [saveRDS] is
#'   called on the entire output as well, saving to `results.rds` by default.
#' * Create an HTML report document summarizing all results.
#'
#' For defaults used in the configuration, see [CFG_DEFAULTS].
#'
#' The workflow above outlines CHIIMP's behavior when called as a standalone
#' program, where [main] loads a configuration file into global options in R and
#' calls [full_analysis].  The public functions linked above can also be used
#' independently; see the documentation and code examples for the individual
#' functions for more information.
#'
#' **The Package structure of the source files, grouped by topic:**
#'  * Main Interface:
#'    * `chiimp.R`: Main entry point for command-line usage ([main]) and R usage
#'    ([full_analysis]).
#'  * Data Analysis:
#'    * `analyze_dataset.R`: High-level interface to analyze all samples
#'      across a given dataset ([analyze_dataset]); used by [full_analysis] to
#'      manage the main part of the processing.
#'    * `summarize_dataset.R`: High-level interface to provide inter-sample
#'      and inter-locus analyses ([summarize_dataset]); used by [full_analysis]
#'      to manage the second stage of the processing.
#'    * `analyze_seqs.R`: Low-level interface to convert raw sequence input
#'      to a data frame of unique sequences ([analyze_seqs]); used
#'      by [analyze_dataset].
#'    * `analyze_sample.R`: Low-level interface to extract per-locus
#'      details from a data frame of unique sequences
#'      ([analyze_sample]); used by [analyze_dataset].
#'    * `summarize_sample.R`: Low-level interface to condense each sample
#'      data frame into a a concise list of consistent attributes, suitable for
#'      binding together across samples for a dataset ([summarize_sample]); used
#'      by [analyze_dataset].
#'    * `categorize.R`: Low-level helper functions used by [summarize_dataset]
#'      for samples with known identity.
#'  * Plotting and reporting:
#'    * `report.R`: Various plotting and summarizing functions used when
#'      rendering a report in [full_analysis].
#'    * `histogram.R`: Sequence histogram plotting tools [histogram]) as used
#'      during [full_analysis].
#'    * `markdown.R`: Various helper functions for adding tables and plots
#'      to an R Markdown report as used in [full_analysis].
#'  * Utility Functions and Configuration:
#'    * `configuration.R`: Configuration handling helper code and the default
#'    configuration options [CFG_DEFAULTS]) used by many chiimp functions.
#'    * `io.R`: various helper input/output functions used loading and
#'      saving sequence data files, spreadsheets, and plots used in multiple
#'      parts of the package.
#'    * `util.R`: Various helper functions used in multiple parts of the
#'      package.
#'
#' @md
#'
"_PACKAGE"

# Analysis ----------------------------------------------------------------

#' Perform a full microsatellite analysis
#'
#' Run all aspects of a microsatellite analysis, using the currently-defined
#' global options (see [load_config] and [apply_config]), and save the
#' corresponding output files.
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @examples
#' \dontrun{
#' # Set up a temporary copy of the CHIIMP test data
#' example_dir <- tempfile()
#' dir.create(example_dir)
#' setwd(example_dir)
#' test_data$write_seqs(test_data$seqs,
#'                      "str-dataset",
#'                       "Replicate1-Sample%s-%s.fasta")
#' locus_attrs_path <- system.file("example_locus_attrs.csv",
#'                                 package = "chiimp")
#' file.copy(locus_attrs_path, "locus_attrs.csv")
#' # Run the example analysis
#' config_path <- system.file("example_config.csv", package = "chiimp")
#' config <- load_config(config_path)
#' apply_config(config)
#' results <- full_analysis()
#' }
#'
#' @export
#' @md
full_analysis <- function() {
  # Only show identifications if a known genotypes table was supplied
  cfg("report_section_identifications", ! is_blank(cfg("genotypes_known")) &&
    cfg("report_section_identifications"))

  logmsg(paste0(
    "Loading dataset: ",
    ifelse(is_blank(cfg("dataset")), cfg("prep_dataset_path"), cfg("dataset")),
    "..."))
  dataset <- if (is_blank(cfg("dataset"))) {
    prepare_dataset()
  } else {
    load_dataset(cfg("fp_dataset"))
  }

  logmsg(paste0("Loading locus attrs: ", cfg("locus_attrs"), "..."))
  locus_attrs <- load_locus_attrs(cfg("locus_attrs"))

  logmsg("Analyzing samples...")
  idx <- match(cfg("sample_summary_func"), sample_summary_funcs, nomatch = 1)
  sample_summary_func <- get(sample_summary_funcs[idx])
  idx <- match(cfg("sample_analysis_func"), sample_analysis_funcs, nomatch = 1)
  sample_analysis_func <- get(sample_analysis_funcs[idx])
  allele_names <- NULL
  if (! is_blank(cfg("allele_names")))
    allele_names <- load_allele_names(cfg("allele_names"))
  results <- analyze_dataset(
    dataset, locus_attrs,
    analysis_opts = list(fraction.min = cfg("min_allele_abundance")),
    summary_opts = list(counts.min = cfg("min_locus_reads")),
    analysis_function = sample_analysis_func,
    summary_function = sample_summary_func,
    known_alleles = allele_names)
  results$allele.names <- allele_names
  results$locus_attrs <- locus_attrs
  logmsg("Summarizing results...")
  genotypes_known <- NULL
  if (! is_blank(cfg("genotypes_known")))
    genotypes_known <- load_genotypes(cfg("genotypes_known"))
  results <- summarize_dataset(results, genotypes_known)
  if (! is_blank(cfg("genotypes_known")))
    results$closest_matches <- find_closest_matches(
      results$dist_mat_known,
      range = cfg("id_dist_range"),
      maximum = cfg("id_dist_max"))
  results$config <- cfg()
  logmsg("Saving output files...")
  save_data(results, results$config)
  if (cfg("report")) {
    logmsg("Creating report...")
    render_report(results, results$config)
  }
  logmsg("Done.")
  return(results)
}

#' Handle full microsatellite analysis from command-line
#'
#' A small wrapper function to read a configuration file path from command-line
#' arguments, load the configuration data (see [load_config]), and run
#' [full_analysis].
#'
#' @param args optional character vector of arguments to use rather than those
#'   detected with [base::commandArgs].
#'
#' @return list of results, with the full configuration list included as
#'   "config."
#'
#' @examples
#' \dontrun{
#' # Set up a temporary copy of the CHIIMP test data
#' example_dir <- tempfile()
#' dir.create(example_dir)
#' setwd(example_dir)
#' test_data$write_seqs(test_data$seqs,
#'                      "str-dataset",
#'                       "Replicate1-Sample%s-%s.fasta")
#' locus_attrs_path <- system.file("example_locus_attrs.csv",
#'                                 package = "chiimp")
#' file.copy(locus_attrs_path, "locus_attrs.csv")
#' # Run the example analysis
#' config_path <- system.file("example_config.csv", package = "chiimp")
#' results <- main(config_path)
#' }
#'
#' @export
#' @md
main <- function(args = NULL) {
  if (missing(args))
    args <- commandArgs(trailingOnly = TRUE)
  desc <- "Identify microsatellite alleles fasta/fastq files"
  p <- argparser::arg_parser(desc)
  p <- argparser::add_argument(p, "config", help = "configuration file path")
  args_parsed <- argparser::parse_args(p, args)
  config_table <- load_config(args_parsed$config)
  apply_config(config_table)
  setwd(dirname(args_parsed$config))
  full_analysis()
}

# Util --------------------------------------------------------------------

#' Render Microsatellite Report Document
#'
#' Write an HTML report summarizing the results of a full microsatellite
#' analysis.
#'
#' @param results list of microsatellite analysis results as produced by
#'   [full_analysis].
#' @param config list of parsed configuration options (see [parse_config]).
#' @md
render_report <- function(results, config) {
  with(config, {
    fp_report_in <- system.file("report", "report.Rmd", package = "chiimp")
    fp_report_out <- file.path(output_path, output_report)
    if (!dir.exists(dirname(fp_report_out)))
      dir.create(dirname(fp_report_out), recursive = TRUE)
    pandoc_metadata <- c(title = report_title,
                         author = report_author,
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
#'   [full_analysis].
#' @param config list of parsed configuration options (see [parse_config]).
#' @md
save_data <- function(results, config) {
  p <- function(thing) file.path(config$output_path, config[[thing]])
  save_histograms(results, p("output_path_histograms"))
  save_results_summary(results$summary, p("output_summary"))
  save_alignments(results$alignments, p("output_path_alignments"))
  save_alignment_images(results$alignments, p("output_path_alignment_images"))
  save_seqfile_data(results$files, p("output_path_processed_files"))
  save_sample_data(results$samples, p("output_path_processed_samples"))
  save_allele_seqs(results$summary, p("output_path_allele_seqs"))
  save_dist_mat(results$dist_mat, p("output_dist_mat"))
  if (! is_blank(config$output_rds))
    saveRDS(results, p("output_rds"))
}

#' Create Pandoc Metadata Argument Strings
#'
#' Convert named list of pandoc metadata options into command-line
#' `--metadata` argument strings.
#'
#' @param metadata named list of
#'   [pandoc metadata options](http://pandoc.org/MANUAL.html#metadata-blocks).
#'
#' @return list of `--metadata` argument strings
#' @md
format_pandoc_args <- function(metadata) {
  metadata <- paste(names(metadata), metadata, sep = ":")
  paste("--metadata=", metadata, sep = "")
}

.getenv <- function() {
  parent.env(parent.frame())
}

#' Set up CHIIMP package environment
#'
#' This loads [CFG_DEFAULTS] from a CSV file when the chiimp package is loaded,
#' and applies the configuration options globally.
#' @md
setup_package <- function() {
  pkg_env <- .getenv()
  cfg_table <- load_config_csv(
    system.file("extdata", "config_defaults.csv", package = "chiimp"))
  for (item in c("CFG_DEFAULTS", "test_data")) {
    if (exists(item, pkg_env)) {
      unlockBinding(item, pkg_env)
    }
  }
  assign("CFG_DEFAULTS", cfg_table, pos = pkg_env)
  apply_config(CFG_DEFAULTS)
  if (! exists("test_data")) {
    assign("test_data", make_helper_data(), pos = pkg_env)
  }
}

.onLoad <- function(libname, pkgname) {
  setup_package()
}
