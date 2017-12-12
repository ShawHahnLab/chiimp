# test load_locus_attrs ---------------------------------------------------

# Take a version of the raw locus_attrs text from helper_data.R, save it to a
# temporary file in TSV format, and return the path.
write_locus_attrs <- function(txt) {
  fp <- tempfile()
  write(gsub(' +', '\t', txt), file = fp)
  fp
}

# create directory of empty files following a given name pattern.
setup_data_dir <- function(replicates, samples, loci, ord=c(1, 2, 3)) {
  dp <- tempfile()
  dir.create(dp)
  f <- function(...) paste(..., sep='-')
  args <- list(replicates, samples, loci)
  args <- args[ord]
  stems <- apply(do.call(expand.grid, args),
                 1,
                 paste,
                 collapse="-")
  fps <- file.path(dp, paste0(stems, '.fastq'))
  lapply(fps, function(fp) cat('', file=fp))
  return(list(dp=dp, fps=fps))
}

test_that("load_locus_attrs parses locus details", {
  # Write the raw text to a temp file, read it back in, and verify that the data
  # structure matches the predefined one.
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
  txt.wrong <- gsub('LengthMin', 'length_min', txt.locus_attrs)
  fp <- write_locus_attrs(txt.wrong)
  expect_warning({
    locus_attrs_test <- load_locus_attrs(fp)
  })
  file.remove(fp)
  expect_equal(nrow(locus_attrs_test), nrow(locus_attrs))
  expect_equal(ncol(locus_attrs_test), ncol(locus_attrs))
  expect_equal(locus_attrs_cols[-1], colnames(locus_attrs_test)[-1])
  expect_equal(c("A", "B", "1", "2"), rownames(locus_attrs_test))
  expect_true(all(locus_attrs == locus_attrs_test))
})

# test prepare_dataset ----------------------------------------------------

test_that("prepare_dataset parses file paths", {
  ## Basic test of prepare_dataset
  replicates <- 1:3
  samples <- 1:5
  loci <- sort(rownames(locus_attrs))
  # by default the field ordering is assumed to be replicate, sample, locus
  data <- setup_data_dir(replicates, samples, loci)
  dataset <- prepare_dataset(data$dp, '(\\d+)-(\\d+)-([A-Za-z0-9]+)')
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
  dataset <- prepare_dataset(data$dp, '([A-Za-z0-9]+)-(\\d+)-(\\d+)', ord)
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
  # Whoops, we left out the loci field in the pattern.  We should get a warning.
  data <- setup_data_dir(replicates, samples, loci)
  expect_warning(dataset <- prepare_dataset(data$dp, '(\\d+)-(\\d+)'))
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
    dataset <- prepare_dataset(data$dp, "(\\d+)-(\\d+)-([A-Za-z0-9]+)")
    }, "Some replicate/sample/locus combinations match multiple files")
})

test_that("prepare_dataset works on nested directories", {
  skip("test not yet implemented")
})
