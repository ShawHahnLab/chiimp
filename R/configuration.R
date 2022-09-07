#' Default microsatellite configuration
#'
#' The entries in this list show the build-time configuration defaults for all
#' aspects of the microsatellite analysis.  These can be overridden by passing a
#' list to \code{\link{full_analysis}} with entries of the same names, or via
#' the configuration file passed to \code{\link{main}} and read via
#' \code{\link{load_config}}.  Check the contents of \code{config.defaults}
#' itself to see all of the build-time defaults.
#'
#' Notable Options:
#'   * \code{dataset_opts}:
#'     * \code{dp}: directory path to input sequence files
#'     * \code{pattern}: regular expression for the input filename pattern
#'     * \code{ord}: order of fields Replicate, Sample, and Locus in in the
#'     input filename pattern.  For example, if Locus is the first field
#'     followed by Replicate and Sample, set \code{ord=c(3, 1, 2)}.
#'   * \code{output}:
#'     * \code{dp}: directory path for saving output data
#'   * \code{fp_dataset}: file path to table of sample attributes to use, rather
#'     than detecting via dataset_opts
#'   * \code{fp_locus_attrs}: file path to locus attributes CSV file
#'   * \code{fp_genotypes_known}: file path to known genotypes CSV file
#' @md
#'
#' @export
config.defaults <- list(
  ## Input and output paths
  fp_dataset = NULL,
  fp_locus_attrs = "locus_attrs.csv",
  fp_allele_names = NULL,
  fp_genotypes_known = NULL,
  ## Arguments to prepare_dataset
  dataset_opts = list(
    dp = "str-data",
    pattern = "(\\d+)-(\\d+)-([A-Za-z0-9]+).fast[aq](?:\\.gz)",
    ord = c(1, 2, 3),
    autorep = FALSE
  ),
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
    dp_processed_files = "processed-files",  # Per-file data tables
    dp_processed_samples = "processed-samples",  # Sample data tables
    dp_allele_seqs = "allele-sequences"  # FASTA sequences for alleles
  ),
  ## Options for analyze_dataset
  dataset_analysis = list(
    # Number of CPU cores to use while processing samples in parallel.  Set to 0
    # to automatically select a number, or 1 to disable parallel execution.
    ncores = 0,
    # Length of suffix on automated allele names in tables.
    name_args = list(hash_len = 6)
  ),
  ## Sample genotyping settings
  seq_analysis = list(
    nrepeats = 3,
    stutter.count.ratio_max = 1 / 3,
    artifact.count.ratio_max = 1 / 3,
    use_reverse_primers = FALSE,
    reverse_primer_r1 = TRUE
  ),
  sample_analysis_func = "analyze_sample",
  sample_analysis_opts = list(fraction.min = 0.05),
  sample_summary_func = "summarize_sample",
  sample_summary_opts = list(counts.min = 500),
  ## Report generation settings
  # Should a report be generated?
  report = TRUE,
  # Should the code executed in the report generation be included in the report?
  report.echo = FALSE,
  # Title and other metadata at top of the report.
  report.title = "Microsatellite Report",
  report.author = NULL,
  # List of vectors of locus names to use to break up tables into reasonable
  # sizes and/or order locus names explicitly.
  report.locus_chunks = NULL,
  # Group together genotype table rows by sample?
  report.group_samples = FALSE,
  # Text to use for NA entries in Replicates column of tables (i.e. Pooled)
  report.na.replicates = "",
  # Text to use for NA allele enries in tables (when the Sample+Locus
  # combination was not present at all, as opposed to when a genotype could not
  # be determined)
  report.na.alleles = "",
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
