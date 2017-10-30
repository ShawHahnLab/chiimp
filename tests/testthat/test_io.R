# test load_locus_attrs ---------------------------------------------------

# Take a version of the raw locus_attrs text from helper_data.R, save it to a
# temporary file in TSV format, and return the path.
write_locus_attrs <- function(txt) {
  fp <- tempfile()
  write(gsub(' +', '\t', txt), file = fp)
  fp
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
  skip("test not yet implemented")
})

test_that("prepare_dataset works on nested directories", {
  skip("test not yet implemented")
})

test_that("prepare_dataset handles different field ordering", {
  skip("test not yet implemented")
})

test_that("prepare_dataset handles broken patterns", {
  skip("test not yet implemented")
})
