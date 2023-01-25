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
    reverse_primer_r1 = TRUE,
    max_mismatches = 0,
    primer_action = "none",
    # individual forward and reverse settings will default to the above values
    # if left as NULL
    max_mismatches_fwd = NULL,
    max_mismatches_rev = NULL,
    primer_action_fwd = NULL,
    primer_action_rev = NULL
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


# NEW ---------------------------------------------------------------------

# spec for any config table
CONFIG_TABLE_SPEC <- data.frame(
  Column = c(
    "Key", "Value", "Description", "Example", "Parser", "OldName", "Comments"),
  Class = "character",
  Required = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)

apply_spec <- function(obj, spec) {
  # TODO check obj against spec
  obj
}

# parse a single character string as a vector of integers separated by
# non-digits
as_integer_vec <- function(txt) {
  if (length(txt) != 1) {
    stop(paste("txt should be of length 1; received", length(txt)))
  }
  as.integer(strsplit(txt, "[^-0-9]+")[[1]])
}

# parse a single character string as a named list of vectors
# name=item1/item2/item3;name2=item4/item5;...
as_locus_vecs <- function(txt) {
  if (length(txt) != 1) {
    stop(paste("txt should be of length 1; received", length(txt)))
  }
  chunks <- strsplit(txt, "; *")[[1]]
  chunk_names <- sub("=.*", "", chunks)
  vecs <- lapply(chunks, function(chunk) {
    strsplit(sub(".*=", "", chunk), "/")[[1]]
  })
  names(vecs) <- chunk_names
  vecs
}

# return ">", "<", or "=" to signify if the first version string given is later,
# earlier, or equivalent to the second version string.
# e.g. "1.2.0" > "1.1.5" 
version_compare <- function(ver1txt, ver2txt) {
  ver1 <- strsplit(ver1txt, "\\.")[[1]]
  ver2 <- strsplit(ver2txt, "\\.")[[1]]
  len <- max(length(ver1), length(ver2))
  ver_mat <- matrix(as.integer(c(ver1[1:len], ver2[1:len])), ncol = 2)
  ver_mat[is.na(ver_mat)] <- 0
  for (idx in seq_len(nrow(ver_mat))) {
    if (ver_mat[idx, 1] > ver_mat[idx, 2]) {
      return(">")
    }
    if (ver_mat[idx, 1] < ver_mat[idx, 2]) {
      return("<")
    }
  }
  return("=")
}

# Does the given vector look "blank"?
# NULL: TRUE (caught by length of 0)
# NA: TRUE
# "": TRUE
# any combo of NA and "": TRUE
# all else: FALSE
is_blank <- function(vec) {
  length(vec) == 0 || all(vec %in% c(NA, ""))
}

# Load a configuration table from CSV into a data frame
# (for now it's still just a dumb grid of character data, not parsed structured
# data)
load_config_table <- function(fp) {
  cfg_table <- read.csv(
    fp, row.names = NULL, colClasses = "character", stringsAsFactors = FALSE)
  cfg_table <- apply_spec(cfg_table, CONFIG_TABLE_SPEC)
  if (exists("CFG_DEFAULTS")) {
    # If we have the defaults loaded (if not, this should *be* loading the
    # defaults) check this file against those values
    config_check_keys(cfg_table)
    config_check_version(cfg_table)
  }
  cfg_table
}

config_check_keys <- function(cfg_table) {
  unknown_txt <- cfg_table$Key[! cfg_table$Key %in% CFG_DEFAULTS$Key]
  if (length(unknown_txt) > 0) {
    warning(paste0(
      "unrecognized config file entries:\n",
      paste(gsub("^", "  ", unknown_txt), collapse = "\n")))
  }
}

config_check_version <- function(cfg_table) {
  version <- cfg_table$Value[match("version", cfg_table$Key)]
  version_pkg <- CFG_DEFAULTS$Value[match("version", CFG_DEFAULTS$Key)]
  ver_cmp <- version_compare(version, version_pkg)
  if (ver_cmp == ">") {
    warning(paste0(
      "Config file version (", version,
      ") > package config version (", version_pkg, ")"))
  }
}

# take a config table, parse each setting, and assign it to a CFG environment.
# take a config table, parse each setting, and return a list of key/value pairs.
parse_config <- function(cfg_table) {
  # if (is.null(cfg_env)) {
  #   cfg_env <- new.env()
  # }
  cfg_list <- list()
  for (idx in seq_len(nrow(cfg_table))) {
    key <- cfg_table$Key[idx]
    val <- cfg_table$Value[idx]
    idx_default <- match(key, CFG_DEFAULTS$Key)
    funcname <- cfg_table$Parser[idx]
    if (is.na(idx_default)) {
      warning(paste("unrecognized config file entry:", key))
    }
    if (is.null(funcname)) {
      # If this config didn't supply a parser, use the default's
      funcname <- CFG_DEFAULTS$Parser[idx_default]
    }
    if (is.na(funcname)) {
      # if the default didn't supply a parser (should only happen for
      # unrecognized entries) use as.character.
      funcname <- "as.character"
    }
    func <- get(funcname)
    cfg_list[[key]] <- func(val)
  }
  cfg_list
}

# just a helper to get from cfg("thing") to options("chiimp.thing"),
# to get a chiimp option, or to get a list of all chiimp options
cfg <- function(key=NULL, val=NULL) {
  if (is.null(key)) {
    opts <- options()
    opts <- opts[grepl("^chiimp\\.", opts)]
    names(opts) <- sub("^chiimp\\.", "", names(opts))
    return(opts)
  }
  if (is.null(val)) {
    return(getOption(paste0("chiimp.", key)))
  }
  cfg_list <- list(val)
  names(cfg_list) <- paste0("chiimp.", key)
  do.call(options, cfg_list)
  val
}

# take a parsed config list and put each entry in options(), with a "chiimp."
# prefix
apply_config <- function(cfg_list, keep = TRUE) {
  if (! keep) {
    options()
    opts <- options()
    opts <- names(opts)[grepl("^chiimp\\.", names(opts))]
    names(opts) <- opts
    opts <- lapply(opts, function(item) NULL)
    do.call(options, opts)
  }
  names(cfg_list) <- paste0("chiimp.", names(cfg_list))
  do.call(options, cfg_list)
}

#' Default Configuration
#' 
#' The default configuration
#' 
#' @name CFG_DEFAULTS
#' @export CFG_DEFAULTS
NULL

#' CHIIMP Configuration
#' 
#' The live configuration.
#' 
#' @name CFG
#' @export CFG
NULL


.getenv <- function() {
  parent.env(parent.frame())
}

setup_package <- function() {
  pkg_env <- .getenv()
  cfg_table <- load_config_table(
    system.file("extdata", "config_defaults.csv", package = "chiimp"))
  cfg_env <- NULL
  for (item in c("CFG_DEFAULTS", "test_data")) {
    if (exists(item, pkg_env)) {
      unlockBinding(item, pkg_env)
    }
  }
  assign("CFG_DEFAULTS", cfg_table, pos = pkg_env)
  cfg_list <- parse_config(CFG_DEFAULTS)
  apply_config(cfg_list)
  if (! exists("test_data")) {
    assign("test_data", make_helper_data(), pos = pkg_env)
  }
}

.onLoad <- function(libname, pkgname) {
  setup_package()
}