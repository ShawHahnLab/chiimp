# Save and load data.

# Load Inputs -------------------------------------------------------------

# Expected column names for locus_attrs
locus_attrs_cols <- c("Locus", "LengthMin", "LengthMax", "LengthBuffer",
                      "Motif", "Primer", "ReversePrimer")

# Exepcted column names for genotypes.known
genotypes_cols <- c("Name", "Locus", "Allele1Seq", "Allele2Seq")

# Expected column names for allele.names
allele_names_cols <- c("Locus", "Name", "Seq")

# Expected column names for dataset
dataset_cols <- c("Filename", "Replicate", "Sample", "Locus")

#' Load configuration file
#'
#' Load a YAML-formatted text file of configuration options for microsatellite
#' analysis.  This is currently just a wrapper around
#' \code{\link[yaml]{yaml.load}}.  The \code{\link{main}} function loads
#' configuration options with this for use by \code{\link{full_analysis}}.
#'
#' @param fp path to configuration file.
#'
#' @return list of configuration options
#'
#' @export
#'
#' @examples
#' filename <- system.file("example_config.yml", package = "chiimp")
#' config <- load_config(filename)
#' # And then: full_analysis(config)
load_config <- function(fp) {
  if (is.na(fp))
    return(list())
  text <- readChar(fp, file.info(fp)$size)
  yaml::yaml.load(text)
}

#' Load and save tables from CSV
#'
#' Load/save a comma-separated table from/to a CSV file.  (These are generic
#' wrapper functions used by more specific loaders like
#' \code{\link{load_locus_attrs}}.)
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to \code{\link[utils]{read.table}} or
#'   \code{\link[utils]{write.table}}.
#'
#' @return data frame
#' 
#' @describeIn load_csv Load CSV
#'
#' @export
load_csv <- function(fp, ...) {
  data <- utils::read.table(fp,
                            header = TRUE,
                            sep = ",",
                            stringsAsFactors = FALSE,
                            ...)
  rownames(data) <- make_rownames(data)
  data
}

#' @describeIn load_csv Save CSV
#' @param data data frame to save to CSV file
#' @export
save_csv <- function(data, fp, ...) {
  mkparents(fp)
  utils::write.table(data,
                     file = fp,
                     sep = ",",
                     na = "",
                     row.names = FALSE,
                     ...)
  data
}

#' Load table of locus attributes
#'
#' Load a comma-separated table of locus attributes to use for analysis.  This
#' is called automatically during \code{\link{full_analysis}}, with the data
#' frame then used by \code{\link{analyze_seqs}} within
#' \code{\link{analyze_dataset}}.
#'
#' @details
#' Columns Required:
#'   * \code{Locus}: Unique identifier for a given locus
#'   * \code{LengthMin}: Minimum known allele sequence length for this locus
#'   * \code{LengthMax}: Minimum known allele sequence length for this locus
#'   * \code{LengthBuffer}: Additional sequence length below \code{LengthMin}
#'     and above \code{LengthMax} to accept for a candidate allele
#'   * \code{Primer}: The forward PCR primer sequence for a given locus, used
#'   when matching sequences to loci
#'   * \code{ReversePrimer}: The reverse PCR primer sequence
#' @md
#'
#' @param fp.locus_attrs path to text file.
#' @param ... additional arguments passed to \code{\link{load_csv}}.
#'
#' @return data frame of locus attributes
#'
#' @export
#'
#' @examples
#' filename <- system.file("example_locus_attrs.csv", package = "chiimp")
#' locus_attrs <- load_locus_attrs(filename)
load_locus_attrs <- function(fp.locus_attrs, ...) {
  data <- load_csv(fp.locus_attrs, ...)
  col.missing <- is.na(match(locus_attrs_cols, colnames(data)))
  if (any(col.missing)) {
    warning(paste("Missing columns in locus_attrs table:",
                  locus_attrs_cols[col.missing]))
  }
  rows_repeated <- data$Locus[-match(data$Locus, data$Locus)]
  if (length(rows_repeated) > 0) {
    warning(paste("Repeated loci in attrs table:", rows_repeated))
  }
  rownames(data) <- make_rownames(data)
  data
}

#' Load table of known allele names
#'
#' Load a comma-separated table of allele names for use in reporting.
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to \code{\link{load_csv}}.
#'
#' @return data frame of allele names
#'
#' @export
load_allele_names <- function(fp, ...) {
  data <- load_csv(fp, ...)
  col.missing <- is.na(match(allele_names_cols, colnames(data)))
  if (any(col.missing)) {
    warning(paste("Missing columns in allele.names table:",
                  allele_names_cols[col.missing]))
  }
  data
}

#' Load table of genotypes
#'
#' Load a comma-separated table of genotypes, one pair of alleles per row.  This
#' information is used to compare samples to genotypes of known individuals in
#' \code{\link{summarize_dataset}}.
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to \code{\link{load_csv}}.
#'
#' @return data frame of genotypes
#'
#' @export
load_genotypes <- function(fp, ...) {
  data <- load_csv(fp, colClasses = "character", na.strings = "", ...)
  col.missing <- is.na(match(genotypes_cols, colnames(data)))
  if (any(col.missing)) {
    warning(paste("Missing columns in genotypes table:",
                  genotypes_cols[col.missing]))
  }
  data <- data[order_entries(data), ]
  rownames(data) <- make_rownames(data)
  data
}

#' Load table of sample attributes
#'
#' Load a comma-separated table of sample attributes for the dataset to be
#' analyzed.  Alternatively, use \code{\link{prepare_dataset}} to automatically
#' read sample attributes from filenames.  If more than one locus is to be
#' analyzed from a single sequencer sample (i.e., multiplexed samples), either
#' the \code{locusmap} argument to \code{prepare_dataset} can be used, or
#' \code{load_dataset} with an explicit mapping of loci to files.
#'
#' @details
#' Columns Required:
#' * Filename: path to each data file
#' * Replicate: identifier for repeated samples; use blanks if not applicable
#' * Sample: identifier for a given biological sample
#' * Locus: locus identifier matching that used in the locus attributes table
#'   (see \code{\link{load_locus_attrs}})
#' @md
#'
#' @param fp path to text file.
#' @param ... additional arguments passed to \code{\link{load_csv}}.
#'
#' @return data frame of sample attributes for the dataset.
#'
#' @export
load_dataset <- function(fp, ...) {
  data <- load_csv(fp, colClasses = "character", na.strings = "", ...)
  col.missing <- is.na(match(dataset_cols, colnames(data)))
  files.missing <- ! file.exists(data$Filename)
  if (any(files.missing)) {
    logmsg(paste("WARNING: Missing", sum(files.missing), "of",
                 length(files.missing), "data files"))
  }
  if (any(col.missing)) {
    warning(paste("Missing columns in dataset table:",
                  dataset_cols[col.missing]))
  }
  rownames(data) <- make_rownames(data)
  data
}

#' Save table of sample attributes
#'
#' Save a comma-separated table of sample attributes.  (This is a convenience
#' function not used automatically in the analysis.)
#'
#' @param data data frame of sample attributes as produced by
#'   \code{\link{prepare_dataset}} or \code{\link{load_dataset}}.
#' @param fp path to text file.
#' @param ... additional arguments passed to \code{\link{save_csv}}.
#'
#' @export
save_dataset <- function(data, fp, ...) {
  save_csv(data, fp, ...)
}

#' Extract Sample Attributes from Filenames
#'
#' Find files matching a pattern in a given directory, and build a data frame
#' of standard sample attributes from fields in the filenames.  Nested
#' directory structures are supported.  Alternatively, use
#' \code{\link{load_dataset}} to load a spreadsheet of sample attributes
#' explicitly.  \code{load_dataset} can be used for cases where more than one
#' locus is to be analyzed from a single sequencer sample (i.e., multiplexed
#' samples), though the \code{locusmap} argument here can allow automatic
#' matching of locus names for multiplexed samples.  If the directory path
#' given does not exist or if no matching files are found, an error is thrown.
#'
#' @param dp directory path to search for matching data files.
#' @param pattern regular expression to use for parsing filenames.  There should
#'   be exactly three groups in the pattern, for Replicate, Sample, and Locus.
#' @param ord integer vector giving order of the fields Replicate, Sample, and
#'   Locus in filenames.  For example, if Locus is the first field followed by
#'   Replicate and Sample, set \code{ord=c(3, 1, 2)}.
#' @param autorep logical allowing for automatic handling of any duplicates
#'   found, labeling them as replicates.  FALSE by default.
#' @param locusmap list of character vectors, each list item name being the
#'   locus text given in the filenames, and each vector being a set of separate
#'   locus names.  Each entry with a locus name text matching one of these list
#'   items will be replaced in the final output with several separate entries,
#'   one for each locus name in the corresponding vector.  (For example,
#'   \code{locusmap=list(ABCD=c("A", "B", "C", "D"))} would take a filename with
#'   "ABCD" in the locus field and split it out into four entries for the four
#'   loci.)
#'
#' @return data frame of metadata for all files found
#'
#' @export
prepare_dataset <- function(dp, pattern, ord = c(1, 2, 3), autorep=FALSE,
                            locusmap=NULL) {
  if (! dir.exists(dp)) {
    stop(paste("ERROR: directory path for data files does not exist:",
               dp))
  }
  # get all matching filenames and extract substrings
  seq_files <- list.files(path = dp,
                          pattern = pattern,
                          full.names = TRUE,
                          recursive = TRUE,
                          include.dirs = FALSE)
  if (! length(seq_files)) {
    stop(paste("ERROR: no data files found:",
               dp))
  }
  seq_file_attrs <- stringr::str_match_all(seq_files, pattern)
  if (! all(sapply(seq_file_attrs, length) == length(ord) + 1)) {
    warning("Some filenames did not match the given pattern")
  }
  # transpose so the list is grouped by filename/attribute rather than
  # entry by entry.
  seq_file_attrs <- lapply(seq(length(ord) + 1), function(y) {
      sapply(seq_file_attrs, "[", y)
    })
  # filename as the first column, then whatever order the attrs are in
  n <- c("Filename", "Replicate", "Sample", "Locus")
  names(seq_file_attrs) <- n[c(1, 1 + ord)]
  data <- do.call(data.frame, c(seq_file_attrs, stringsAsFactors = FALSE))
  data$Filename <- seq_files
  data$Replicate <- as.integer(ifelse(data$Replicate == "",
                                      NA,
                                      data$Replicate))

  # If specified, map the locus text into multiple loci per sample (e.g.
  # multiplexed)
  if (! is.null(locusmap)) {
    data <- do.call(rbind, lapply(1:nrow(data), function(i) {
      col_locus <- match("Locus", colnames(data))
      cols <- as.list(data[i, -col_locus])
      locus <- data[i, col_locus]
      locus_new <- locusmap[[locus]]
      if (is.null(locus_new)) {
        locus_new <- locus
      }
      args <- c(cols,
              Locus = list(locus_new),
              stringsAsFactors = FALSE)
      do.call(data.frame, args)
    }))
  }

  # order by locus/sample/replicate
  data <- data[order_entries(data), ]

  # If specified, automatically number duplicates as replicates, if they don't
  # have a replicate number already.
  if (autorep) {
    data <- do.call(rbind, lapply(split(data,
                                           paste(data$Sample,
                                                 data$Locus)),
           function(chunk) {
             chunk$Replicate <- ifelse(is.na(chunk$Replicate),
                                       1:nrow(chunk),
                                       chunk$Replicate)
             chunk
    }))
    data <- data[order_entries(data), ]
  }

  rownames(data) <- make_rownames(data)
  # complain if any replicate/sample/locus combo matches more than one entry
  if (max(table(paste(data$Replicate, data$Sample, data$Locus))) > 1) {
    warning("Some replicate/sample/locus combinations match multiple files")
  }
  return(data)
}

#' Load vector of sequences from FASTA/FASTQ file
#'
#' Load a vector of character sequences from the given path.  This is just a
#' wrapper around \code{\link[dnar:read.fa]{dnar}} to choose the parser based on
#' filename.  Only the sequences are returned, not IDs or quality scores.
#'
#' @param fp path to sequence file
#'
#' @return vector of sequences
#'
#' @export
load_seqs <- function(fp) {
  if (length(grep("q$", fp)) || length(grep("q.gz$", fp)))
    loadfunc <- dnar::read.fastq
  else
    loadfunc <- dnar::read.fa
  loadfunc(fp)$seq
}

# Output Saving -----------------------------------------------------------

#' Save dataset summary to text file
#'
#' Save the dataset summary produced by \code{\link{analyze_dataset}} to the
#' specified file path in CSV format.
#'
#' @param results_summary summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param fp output file path.
#'
#' @export
save_results_summary <- function(results_summary, fp) {
  save_csv(results_summary, fp)
}

#' Save identified alleles to FASTA files
#'
#' Take the alleles identified by \code{\link{analyze_dataset}} in the summary
#' data frame and save each entry to a separate FASTA file.  Samples identified
#' as homozygous will have one sequence written rather than two.  Entries with
#' no identified alleles will be skipped.
#'
#' @param results_summary summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param dp output directory path to use for all files.
#'
#' @export
save_allele_seqs <- function(results_summary, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  results_summary$.row <- rownames(results_summary)
  invisible(apply(results_summary, 1, function(r) {
    seqs <- r[c("Allele1Seq", "Allele2Seq")]
    # Ignore NA entries, and skip this sample if there are no sequences
    seqs <- seqs[!is.na(seqs)]
    if (length(seqs) == 0)
      return()
    # Create sequence names with a unique ID based on the sample and some other
    # stats for convenience
    seq.names <- sapply(seq_along(seqs), function(i) {
      seq.id <- paste(r[".row"], i, sep = "_")
      a.len <- r[paste0("Allele", i, "Length")]
      a.cts <- r[paste0("Allele", i, "Count")]
      seq.extra <- paste(paste0(a.len, "bp"),
                         paste0("allele/locus/total=",
                                paste(as.integer(a.cts),
                                      as.integer(r["CountLocus"]),
                                      as.integer(r["CountTotal"]),
                                      sep = "/")))
      paste(seq.id, seq.extra)
    })
    fp <- file.path(dp, paste0(r[".row"], ".fasta"))
    dnar::write.fa(seq.names, seqs, fp)
  }))
}

#' Save per-file processed data to text files
#'
#' Save each per-file data frame produced by \code{\link{analyze_dataset}} (via
#' \code{\link{analyze_seqs}}) to a separate file in the specified directory
#' path, in CSV format.  The directory structure will start at the first shared
#' directory of the input file paths.
#' For example, if the inputs were \code{/data/run1/file.fastq} and
#' \code{/data/run2/file.fastq} there will be run1 and run2 directories inside
#' the given \code{dp} directory.
#'
#' @param results_file_data list of per-file data frames as produced by
#'   \code{\link{analyze_dataset}}.
#' @param dp output directory path to use for all files.
#'
#' @export
save_seqfile_data <- function(results_file_data, dp) {
  fps_rel <- remove_shared_root_dir(names(results_file_data))
  invisible(lapply(names(results_file_data), function(n) {
    fp_this <- fps_rel[n]
    fp <- file.path(dp, paste0(fp_this, ".csv"))
    save_csv(results_file_data[[n]], fp, quote = FALSE)
  }))
}

#' Save per-sample processed data to text files
#'
#' Save each per-sample data frame produced by \code{\link{analyze_dataset}}
#' (via \code{\link{analyze_sample}}) to a separate file in the specified
#' directory path, in CSV format.
#'
#' @param results_data list of per-sample data frames as produced by
#'   \code{\link{analyze_dataset}}.
#' @param dp output directory path to use for all files.
#'
#' @export
save_sample_data <- function(results_data, dp) {
  invisible(lapply(names(results_data), function(n) {
    fp <- file.path(dp, paste0(n, ".csv"))
    save_csv(results_data[[n]], fp, quote = FALSE)
  }))
}

#' Save alignments to FASTA files
#'
#' Take a list of alignments, one per locus, and save each to a separate FASTA
#' file in a specified directory.  If any of the per-locus alignment objects is
#' NA it will be skipped.  These are produced by \code{\link{summarize_dataset}}
#' via \code{\link{align_alleles}}.
#'
#' @param alignments list of MSA alignment objects, such as created by
#'   \code{\link{summarize_dataset}} via \code{\link{align_alleles}}.  The name
#'   of each alignment will be used for its filename.
#' @param dp output directory path.
#'
#' @export
save_alignments <- function(alignments, dp) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if (!is.null(alignments[[loc]])) {
      dna <- as.character(alignments[[loc]])
      fp <- file.path(dp, paste0(loc, ".fasta"))
      dnar::write.fa(names = names(dna),
                     dna = dna,
                     fileName = fp)
    }
  }))
}

#' Save alignment visualizations to image files
#'
#' Take a list of alignments, one per locus, and save a plot of each to a
#' separate image file in a specified directory.  If any of the per-locus
#' alignment objects is NA it will be skipped.  These are produced by
#' \code{\link{summarize_dataset}} via \code{\link{align_alleles}}.
#'
#' @param alignments list of MSA alignment objects, such as created by
#'   \code{\link{summarize_dataset}} via \code{\link{align_alleles}}.  The name
#'   of each alignment will be used for its filename.
#' @param dp output directory path.
#' @param image.func name of function to call for saving each image.
#' @param width integer width of image.
#' @param height integer height of image.
#' @param res integer resolution of image in PPI.
#'
#' @export
save_alignment_images <- function(alignments, dp, image.func="png",
                                  width=1600, height=1200, res=150) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(alignments), function(loc) {
    if (!is.null(alignments[[loc]])) {
      fp <- file.path(dp, paste(loc, image.func, sep = "."))
      img.call <- call(image.func,
                       fp,
                       width = width,
                       height = height,
                       res = res)
      eval(img.call)
      plot_alignment(alignments[[loc]],
                     main = paste("Alignment for Locus", loc))
      grDevices::dev.off()

    }
  }))
}

#' Save sequence histogram visualizations to image files
#'
#' Take a full results list and save a histogram (via \code{\link{histogram}})
#' of each sample to a separate image file in a specified directory.
#'
#' @param results list of results as created by \code{\link{analyze_dataset}}.
#' @param dp output directory path.
#' @param image.func name of function to call for saving each image.
#' @param width integer width of image.
#' @param height integer height of image.
#' @param res integer resolution of image in PPI.
#'
#' @export
save_histograms <- function(results, dp, image.func="png",
                            width=1600, height=1200, res=150) {
  if (!dir.exists(dp))
    dir.create(dp, recursive = TRUE)
  invisible(lapply(names(results$samples), function(entry) {
    fp <- file.path(dp, paste(entry, image.func, sep = "."))
    img.call <- call(image.func,
                     fp,
                     width = width,
                     height = height,
                     res = res)
    fn <- results$summary[entry, "Filename"]
    seq_data <- results$files[[fn]]
    sample_data <- results$samples[[entry]]
    eval(img.call)
    histogram(seq_data = seq_data,
              sample_data = sample_data,
              main = entry)
    grDevices::dev.off()
  }))
}

#' Save sample distance matrix to text file
#'
#' Save the inter-sample distance matrix produced by
#' \code{\link{summarize_dataset}} to the specified file path in CSV format.
#'
#' @param dist_mat matrix produced by \code{\link{summarize_dataset}} via
#'   \code{\link{make_dist_mat}}.
#' @param fp output file path.
#'
#' @export
save_dist_mat <- function(dist_mat, fp) {
  save_csv(dist_mat, fp)
}


# Helpers -----------------------------------------------------------------


# Create all parent directories as needed
mkparents <- function(fp) {
  if (!dir.exists(dirname(fp)))
    dir.create(dirname(fp), recursive = TRUE)
}
