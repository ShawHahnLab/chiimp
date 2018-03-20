context("Test input/output functions")


# helpers -------------------------------------------------------------------


# Take a version of the raw locus_attrs text from helper_data.R, save it to a
# temporary file in TSV format, and return the path.
write_locus_attrs <- function(txt) {
  fp <- tempfile()
  write(gsub(" +", ",", txt), file = fp)
  fp
}

# create directory of empty files following a given name pattern.
setup_data_dir <- function(replicates, samples, loci, ord=c(1, 2, 3)) {
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
  lapply(fps, function(fp) cat("", file = fp))
  return(list(dp = dp, fps = fps, pattern = pattern))
}

# simple example dataset data frame
setup_dataset <- function(reps=1:3, samps=1:5,
                          loci=sort(rownames(test_data$locus_attrs))) {
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

with(test_data, {


# test load_locus_attrs ---------------------------------------------------


  test_that("load_locus_attrs parses locus details", {
    # Write the raw text to a temp file, read it back in, and verify that the
    # data structure matches the predefined one.
    fp <- write_locus_attrs(txt.locus_attrs)
    locus_attrs_test <- load_locus_attrs(fp)
    file.remove(fp)
    expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
    expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
    expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
    expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
    expect_true(all(locus_attrs == locus_attrs_test))
  })

  test_that("load_locus_attrs handles missing column names", {
    # It should throw a warning if one of the expected columns is missing, and
    # then proceed.  The other columns should still match up.
    txt.wrong <- gsub("LengthMin", "length_min", txt.locus_attrs)
    fp <- write_locus_attrs(txt.wrong)
    expect_warning({
      locus_attrs_test <- load_locus_attrs(fp)
    })
    file.remove(fp)
    expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
    expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
    expect_equal(locus_attrs_cols[-2], colnames(locus_attrs_test)[-2])
    expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
    expect_true(all(locus_attrs == locus_attrs_test))
  })


  test_that("load_locus_attrs complains for repeated loci", {
    # It should throw a warning if any Locus names are repeated.
    txt.wrong <- gsub("\nB ", "\nA ", txt.locus_attrs)
    fp <- write_locus_attrs(txt.wrong)
    expect_warning({
      locus_attrs_test <- load_locus_attrs(fp)
    })
    file.remove(fp)
    # The data frame will still be constructed, but the second A row will have
    # the row name "A.1".
    expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
    expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
    expect_equal(locus_attrs_cols, colnames(locus_attrs_test))
    expect_equal(c("A", "A.1", "1", "2"), rownames(locus_attrs_test))
    expect_true(all(locus_attrs[, -1] == locus_attrs_test[, -1]))
  })


# test load_dataset -------------------------------------------------------


  test_that("load_dataset loads sample attributes", {
    dataset_known <- setup_dataset()
    fp <- tempfile()
    write.csv(dataset_known, file = fp, na = "", row.names = F)
    dataset <- load_dataset(fp)
    expect_identical(dataset, dataset_known)
  })

# test save_dataset -------------------------------------------------------


  test_that("save_dataset saves sample attributes", {
    # We'll just make sure that saving and then re-loading provides what was
    # saved.  load_dataset is tested separately above.
    dataset_known <- setup_dataset()
    fp <- tempfile()
    save_dataset(dataset_known, fp)
    dataset <- load_dataset(fp)
    expect_identical(dataset, dataset_known)
  })

# test prepare_dataset ----------------------------------------------------


  test_that("prepare_dataset parses file paths", {
    ## Basic test of prepare_dataset
    replicates <- 1:3
    samples <- 1:5
    loci <- sort(rownames(locus_attrs))
    # by default the field ordering is assumed to be replicate, sample, locus
    data <- setup_data_dir(replicates, samples, loci)
    dataset <- prepare_dataset(data$dp, data$pattern)
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
    loci <- sort(rownames(locus_attrs))
    # order: 1=replicate, 2=sample, 3=locus
    # --> Locus, Replicate, Sample
    ord <- c(3, 1, 2)
    data <- setup_data_dir(replicates, samples, loci, ord)
    dataset <- prepare_dataset(data$dp, "([A-Za-z0-9]+)-(\\d+)-(\\d+)", ord)
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
    loci <- sort(rownames(locus_attrs))
    # Whoops, we left out the loci field in the pattern.  We should get a
    # warning.
    data <- setup_data_dir(replicates, samples, loci)
    expect_warning(dataset <- prepare_dataset(data$dp, "(\\d+)-(\\d+)"))
  })

  test_that("prepare_dataset warns of repeated identifier rows", {
    # It should throw a warning if there are multiple rows for any
    # Replicate/Sample/Locus combination.
    replicates <- 1:3
    samples <- 1:5
    loci <- sort(rownames(locus_attrs))
    data <- setup_data_dir(replicates, samples, loci)
    cat("", file = paste0(data$fps[3], ".extra"))
    expect_warning({
        dataset <- prepare_dataset(data$dp, data$pattern)
      },
      "Some replicate/sample/locus combinations match multiple files")
  })

  test_that("prepare_dataset can autolabel replicates", {
    # With autorep=TRUE it should automatically number replicates.
    replicates <- 1
    samples <- 1:5
    loci <- sort(rownames(locus_attrs))
    data <- setup_data_dir(replicates, samples, loci)
    cat("", file = paste0(data$fps[3], ".2"))
    cat("", file = paste0(data$fps[3], ".3"))
    dataset <- prepare_dataset(data$dp,
                               pattern = "()1-(\\d+)-([A-Za-z0-9]+)",
                               autorep = TRUE)
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
    skip("test not yet implemented")
  })

  test_that("prepare_dataset can separate multiplexed samples", {
    replicates <- 1
    samples <- 1:5
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
    # Aside from the different filenames, does everything match up?
    dataset_known$Filename <- dataset$Filename
    expect_equal(dataset, dataset_known)
  })

})
