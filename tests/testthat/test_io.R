testpth <- function(fname) normalizePath(test_path("data", "io", fname))
testrds <- function(fname) readRDS(test_path("data", "io", fname))


# test load_config --------------------------------------------------------


test_that("load_config loads config YAML files", {
  config_path <- testpth("config.yml")
  config <- load_config(config_path)
  expect_equal(
    config,
    list(fp_dataset = "samples.csv", output = list(fp_rds = "results.rds"))
  )
})

test_that("load_config handles unexpected entries", {
  config_path <- testpth("config_unrecognized_key.yml")
  # it should warn about whatever config entries are unknown, but still load and
  # return whatever's there
  expect_warning(
    config <- load_config(config_path),
    paste(
      "unrecognized config file entries:",
      "  unrecognized",
      "  dataset_analysis:name_args:unknown", sep = "\n"))
  expect_equal(
    config,
    list(
      fp_dataset = "samples.csv",
      output = list(fp_rds = "results.rds"),
      unrecognized = 10,
      dataset_analysis = list(name_args = list(unknown = 5)))
  )
})


# test load_csv/save_csv --------------------------------------------------


# We'll use the same locus_attrs as below for these tests.

test_that("load_csv parses CSV files", {
  locus_attrs_test <- load_csv(testpth("locus_attrs.csv"))
  locus_attrs <- testrds("locus_attrs.rds")
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})

test_that("load_csv parses CSV files with unknown columns", {
  # For cases where our CSV file has none of the recognized column names, it
  # should still work but just give us generic rownames and a warning.
  data_expected <- testrds("misc.csv.rds")
  expect_warning(
    data <- load_csv(testpth("misc.csv")),
    "no recognized columns for entry id")
  expect_equal(data, data_expected)
})

test_that("load_csv uses row names given", {
  # If we give a row.names argument it should pass that on to read.table and
  # shouldn't automatically generate row names
  data_expected <- testrds("misc.csv.rds")
  rownames(data_expected) <- NULL
  data <- load_csv(testpth("misc.csv"), row.names = NULL)
  expect_equal(data, data_expected)
})

test_that("save_csv saves CSV files", {
  locus_attrs <- testrds("locus_attrs.rds")
  within_tmpdir({
    fp <- tempfile(tmpdir = ".")
    save_csv(locus_attrs, fp)
    locus_attrs_test <- load_locus_attrs(fp)
  })
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})

test_that("save_csv saves CSV files with unknown columns", {
  # save_csv shouldn't care about the columns but we'll make sure here.
  data_expected <- testrds("misc.csv.rds")
  rownames(data_expected) <- c("entry1", "entry2")
  within_tmpdir({
    fp <- tempfile(tmpdir = ".")
    save_csv(data_expected, fp)
    expect_warning(data <- load_csv(fp))
  })
  expect_equal(data, data_expected)
})

test_that("save_csv makes parent dirs when saving", {
  # the parent directories here don't yet exist, so for this to work, save_csv
  # will need to create them
  locus_attrs <- testrds("locus_attrs.rds")
  within_tmpdir({
    fp <- file.path(tempfile(tmpdir = "."), basename(tempfile()))
    save_csv(locus_attrs, fp)
    locus_attrs_test <- load_locus_attrs(fp)
  })
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})


# test load_locus_attrs ---------------------------------------------------


test_that("load_locus_attrs parses locus details", {
  # Write the raw text to a temp file, read it back in, and verify that the data
  # structure matches the predefined one.
  locus_attrs_test <- load_csv(testpth("locus_attrs.csv"))
  locus_attrs <- testrds("locus_attrs.rds")
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})

test_that("load_locus_attrs handles missing column names", {
  # It should throw a warning if one of the expected columns is missing, and
  # then proceed.  The other columns should still match up.
  locus_attrs <- testrds("locus_attrs.rds")
  expect_warning({
    locus_attrs_test <- load_locus_attrs(testpth("locus_attrs_wrongcol.csv"))
  })
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols[-2], colnames(locus_attrs_test)[-2])
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})

test_that("load_locus_attrs complains for repeated loci", {
  # It should throw a warning if any Locus names are repeated.
  locus_attrs <- testrds("locus_attrs_dups.rds")
  expect_warning({
    locus_attrs_test <- load_locus_attrs(testpth("locus_attrs_dups.csv"))
  })
  # The data frame will still be constructed, but the second A row will have the
  # row name "A.1".
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
  expect_equal(c("A", "A.1", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs[, -1] == locus_attrs_test[, -1]))
})


# test load_dataset -------------------------------------------------------


test_that("load_dataset loads sample attributes", {
  dataset_path <- testpth("dataset.csv")
  dataset_known <- testrds("dataset.rds")
  within_tmpdir({
    # Touch the input files so they at least exist
    touch(as.character(read.csv(dataset_path)$Filename))
    expect_silent(dataset <- load_dataset(dataset_path))
  })
  expect_identical(dataset, dataset_known)
})

test_that("load_dataset warns of missing files", {
  # Here we'll skip writing any actual data files, so load_dataset should
  # complain.
  dataset_path <- testpth("dataset.csv")
  dataset_known <- testrds("dataset.rds")
  # expect_message and capture_messages both do NOT catch text send to stderr,
  # though capture.output(..., type = "message") does.
  msg <- capture.output(dataset <- load_dataset(dataset_path), type = "message")
  expect_true(length(grep("WARNING: Missing 60 of 60 data files", msg)) == 1)
  expect_identical(dataset, dataset_known)
})

test_that("load_dataset warns of duplicated sample/replicate/locus combos", {
  dataset_path <- testpth("dataset_dups.csv")
  dataset_known <- testrds("dataset_dups.rds")
  within_tmpdir({
    touch(as.character(read.csv(dataset_path)$Filename))
    expect_warning(
      dataset <- load_dataset(dataset_path),
      paste(
        "Duplicated sample/replicate/locus entries for",
        "1/1/1, 1/1/2, 1/1/A, 1/1/B"))
  })
  expect_identical(dataset, dataset_known)
})


# test save_dataset -------------------------------------------------------


test_that("save_dataset saves sample attributes", {
  # We'll just make sure that saving and then re-loading provides what was
  # saved.  load_dataset is tested separately above.
  dataset_known <- testrds("dataset.rds")
  within_tmpdir({
    touch(dataset_known$Filename)
    fp <- tempfile(tmpdir = ".")
    save_dataset(dataset_known, fp)
    dataset <- load_dataset(fp)
  })
  expect_identical(dataset, dataset_known)
})


# test save_seqfile_data --------------------------------------------------


test_that("save_seqfile_data saves per-file information", {
  # Just make sure the expected files are created based on the existing
  # filenames.  The names should be the existing names with an added csv
  # extension with a flat directory structure.
  locus_attrs <- testrds("locus_attrs.rds")
  seqs <- testrds("seqs.rds")
  within_tmpdir({
    write_seqs(seqs, "data")
    dataset <- prepare_dataset("data", "()(\\d+)-([A-Za-z0-9]+).fasta")
    results <- analyze_dataset(
      dataset, locus_attrs,
      analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500),
      nrepeats = 3,
      ncores = 1)
    dp_out <- file.path("results", "processed_files")
    save_seqfile_data(results$files, dp_out)
    fps_expected <- sort(file.path(
      dp_out,
      paste0(basename(names(results$files)), ".csv")))
    fps_observed <- sort(list.files(
      dp_out, recursive = TRUE, full.names = TRUE))
  })
  expect_equal(fps_observed, fps_expected)
})

test_that("save_seqfile_data works with directory trees", {
  # In this case we should get two sub-directories in the output.
  locus_attrs <- testrds("locus_attrs.rds")
  seqs <- testrds("seqs.rds")
  within_tmpdir({
    write_seqs(seqs, file.path("data", "set1"))
    write_seqs(seqs, file.path("data", "set2"))
    dataset <- prepare_dataset(
      "data", pattern = "()(\\d+)-([A-Za-z0-9]+).fasta", autorep = TRUE)
    results <- analyze_dataset(
      dataset, locus_attrs,
      analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500),
      nrepeats = 3,
      ncores = 1)
    dp_out <- file.path("results", "processed_files")
    save_seqfile_data(results$files, dp_out)
    fps_expected <- sort(file.path(
      dp_out,
      basename(dirname(names(results$files))),
      paste0(basename(names(results$files)), ".csv")))
    fps_observed <- sort(list.files(
      dp_out, recursive = TRUE, full.names = TRUE))
  })
  expect_equal(fps_observed, fps_expected)
})

test_that("save_seqfile_data works with Windows-style paths", {
  # This should behave the same as the above test despite the backslashes
  # substituted in for path separators.
  locus_attrs <- testrds("locus_attrs.rds")
  seqs <- testrds("seqs.rds")
  within_tmpdir({
    write_seqs(seqs, file.path("data", "set1"))
    write_seqs(seqs, file.path("data", "set2"))
    dataset <- prepare_dataset(
      "data", pattern = "()(\\d+)-([A-Za-z0-9]+).fasta", autorep = TRUE)
    results <- analyze_dataset(
      dataset, locus_attrs,
      analysis_opts = list(fraction.min = 0.05),
      summary_opts = list(counts.min = 500),
      nrepeats = 3,
      ncores = 1)
    dp_out <- file.path("results", "processed_files")
    names(results$files) <- gsub("/", "\\\\", names(results$files))
    save_seqfile_data(results$files, dp_out)
    names(results$files) <- gsub("\\\\", "/", names(results$files))
    fps_expected <- sort(file.path(
      dp_out,
      basename(dirname(names(results$files))),
      paste0(basename(names(results$files)), ".csv")))
    fps_observed <- sort(list.files(
      dp_out, recursive = TRUE, full.names = TRUE))
    fps_expected <- normalizePath(fps_expected)
    fps_observed <- normalizePath(fps_observed)
  })
  # Normalize any lingering \ or / inconsistencies, so this test should also
  # pass on Windows itself.
  expect_equal(fps_observed, fps_expected)
})


# test prepare_dataset ----------------------------------------------------


# create directory of empty files following a given name pattern.
setup_data_dir <- function(replicates, samples, loci, ord = c(1, 2, 3)) {
  dp <- tempfile()
  dir.create(dp)
  args <- list(replicates, samples, loci)
  args <- args[ord]
  stems <- apply(do.call(expand.grid, args),
                 1,
                 paste,
                 collapse = "-")
  # (hardcoded to match stems)
  pattern <- "(\\d+)-(\\d+)-([A-Za-z0-9]+)"
  fps <- file.path(dp, paste0(stems, ".fastq"))
  touch(fps)
  return(list(dp = dp, fps = fps, pattern = pattern))
}

# simple example dataset data frame
setup_dataset <- function(reps = 1:3, samps = 1:5,
                          loci = c(1, 2, "A", "B")) {
  ids <- as.character(outer(reps,
                            outer(samps,
                                  loci,
                                  FUN = paste,
                                  sep = "-"),
                            FUN = paste,
                            sep = "-"))
  dataset <- data.frame(
    Filename = paste0(ids, ".fasta"),
    Replicate = as.character(rep(reps, length(samps) * length(loci))),
    Sample = as.character(rep(samps,
                              times = length(loci),
                              each = length(reps))),
    Locus = rep(loci, each = length(samps) * length(reps)),
    stringsAsFactors = FALSE)
  rownames(dataset) <- with(dataset, paste(Sample,
                                           Replicate,
                                           Locus,
                                           sep = "-"))
  dataset
}


test_that("prepare_dataset parses file paths", {
  ## Basic test of prepare_dataset
  replicates <- 1:3
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  # by default the field ordering is assumed to be replicate, sample, locus
  data <- setup_data_dir(replicates, samples, loci)
  dataset <- prepare_dataset(data$dp, data$pattern)
  unlink(x = data$dp, recursive = TRUE)
  expect_equal(colnames(dataset),
               c("Filename", "Replicate", "Sample", "Locus"))
  expect_equal(sort(dataset$Filename), sort(data$fps))
  # The row ordering will vary loci slowest (largest chunks), then samples,
  # then replicates.
  # Loci first, repeated across samples and replicates
  expect_equal(as.character(dataset$Locus),
               rep(loci, each = 15))
  # Samples, repeated across replicates and then cycled across loci
  expect_equal(as.character(dataset$Sample),
               as.character(rep(rep(samples, each = 3), 4)))
  # Replicates, cycled across samples and loci
  expect_equal(dataset$Replicate,
               rep(replicates, 20))
})

test_that("prepare_dataset handles different field ordering", {
  replicates <- 1:3
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  # order: 1=replicate, 2=sample, 3=locus
  # --> Locus, Replicate, Sample
  ord <- c(3, 1, 2)
  data <- setup_data_dir(replicates, samples, loci, ord)
  dataset <- prepare_dataset(data$dp, "([A-Za-z0-9]+)-(\\d+)-(\\d+)", ord)
  unlink(x = data$dp, recursive = TRUE)
  expect_equal(colnames(dataset),
               c("Filename", "Locus", "Replicate", "Sample"))
  expect_equal(sort(dataset$Filename), sort(data$fps))
  # These should be the same as for the previous test; the ordering of rows
  # shouldn't change.
  expect_equal(as.character(dataset$Locus),
               rep(loci, each = 15))
  expect_equal(as.character(dataset$Sample),
               as.character(rep(rep(samples, each = 3), 4)))
  expect_equal(dataset$Replicate,
               rep(replicates, 20))
})

test_that("prepare_dataset handles broken patterns", {
  replicates <- 1:3
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  # Whoops, we left out the loci field in the pattern.  We should get a
  # warning.
  data <- setup_data_dir(replicates, samples, loci)
  expect_warning(dataset <- prepare_dataset(data$dp, "(\\d+)-(\\d+)"))
  unlink(x = data$dp, recursive = TRUE)
})

test_that("prepare_dataset warns of repeated identifier rows", {
  # It should throw a warning if there are multiple rows for any
  # Replicate/Sample/Locus combination.
  replicates <- 1:3
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  data <- setup_data_dir(replicates, samples, loci)
  touch(paste0(data$fps[3], ".extra"))
  expect_warning({
    dataset <- prepare_dataset(data$dp, data$pattern)
  },
  "Some replicate/sample/locus combinations match multiple files")
  unlink(x = data$dp, recursive = TRUE)
})

test_that("prepare_dataset can autolabel replicates", {
  # With autorep=TRUE it should automatically number replicates.
  replicates <- 1
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  data <- setup_data_dir(replicates, samples, loci)
  touch(paste0(data$fps[3], ".2"))
  touch(paste0(data$fps[3], ".3"))
  dataset <- prepare_dataset(data$dp,
                             pattern = "()1-(\\d+)-([A-Za-z0-9]+)",
                             autorep = TRUE)
  unlink(x = data$dp, recursive = TRUE)
  extras <- paste0(data$fps[3], c(".2", ".3"))
  expect_equal(sort(dataset$Filename), sort(c(data$fps, extras)))
  expect_equal(as.character(dataset$Locus),
               c(loci[1], loci[1], rep(loci, each = 5)))
  s <- as.character(rep(samples, 4))
  s <- c(s[1:3], s[3], s[3], s[4:length(s)])
  expect_equal(as.character(dataset$Sample), s)
  r <- rep(1, nrow(dataset))
  r[4:5] <- 2:3
  expect_equal(dataset$Replicate, r)
})

test_that("prepare_dataset works on nested directories", {
  ## Create two separate sets of files in subdirectories
  replicates <- 1:3
  samples1 <- 1:3
  samples2 <- 4:6
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  data1 <- setup_data_dir(replicates, samples1, loci)
  data2 <- setup_data_dir(replicates, samples2, loci)
  # move both sets into a single parent directory
  dp <- tempfile()
  dir.create(dp)
  file.copy(data1$dp, dp, recursive = TRUE)
  file.copy(data2$dp, dp, recursive = TRUE)
  # build dataset from parent directory
  dataset <- prepare_dataset(dp, data1$pattern)
  expect_equal(colnames(dataset),
               c("Filename", "Replicate", "Sample", "Locus"))
  expect_equal(sort(dataset$Filename),
               sort(list.files(dp, recursive = TRUE, full.names = TRUE)))
  unlink(x = c(data1, data2, dp), recursive = TRUE)
})


test_that("prepare_dataset can separate multiplexed samples", {
  replicates <- 1
  samples <- 1:5
  locus_attrs <- testrds("locus_attrs.rds")
  loci <- sort(rownames(locus_attrs))
  # Set up locus name placeholders for multiplexed samples, lumping loci
  # together two at a time (assumes even number)
  locusmap <- sapply(1:2, function(i) {
      loci[seq(i, length(loci), 2)]
    },
    simplify = FALSE)
  names(locusmap) <- sapply(locusmap, paste, collapse = "")
  # Write data using multiplex version of locus names
  data <- setup_data_dir(replicates, samples, names(locusmap))
  # Create the full dataset table we expect to see
  dataset_known <- setup_dataset(reps = replicates,
                                 samp = samples,
                                 loci = loci)
  # prepare_dataset gives integers for replicate.  Should reconcile that up
  # with load_dataset and setup_dataset.
  dataset_known$Replicate <- as.integer(dataset_known$Replicate)
  # Read dataset from disk using the mapping of locus names
  dataset <- prepare_dataset(data$dp, data$pattern, locusmap = locusmap)
  unlink(x = data$dp, recursive = TRUE)
  # Aside from the different filenames, does everything match up?
  dataset_known$Filename <- dataset$Filename
  expect_equal(dataset, dataset_known)
})

test_that("prepare_dataset handles missing data directory", {
  dp <- tempfile()
  expect_error({
      prepare_dataset(dp, "(\\d+)-(\\d+)-([A-Za-z0-9]+)")
    },
    paste("ERROR: directory path for data files does not exist:", dp),
    fixed = TRUE)
})

test_that("prepare_dataset handles no-samples case", {
  dp <- tempfile()
  dir.create(dp)
  expect_error({
      prepare_dataset(dp, "(\\d+)-(\\d+)-([A-Za-z0-9]+)")
    },
    paste("ERROR: no data files found:", dp),
    fixed = TRUE)
  unlink(x = dp, recursive = TRUE)
})