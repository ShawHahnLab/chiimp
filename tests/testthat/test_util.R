test_that("make_entry_id makes identifiers for STR data frames", {
})

test_that("make_rownames makes unique row names for STR data frames", {
})

test_that("order_entries provides ordering vector for STR data frames", {
})

test_that("make_allele_names autogenerates allele names for sequences", {
})

test_that("order_alleles provides ordering vector for allele name pairs", {
})

test_that("name_alleles_in_table adds allele name cols to STR data frames", {
})

test_that("remove_shared_root_dir removes common parent directories of paths", {
})

test_that("logmsg writes messages to stderr if verbose option is set", {
})


# revcmp ------------------------------------------------------------------


test_that("revcmp reverse complements", {
  # basic cases
  expect_identical(revcmp("GCCA"), "TGGC")
  expect_identical(revcmp("GCca"), "tgGC")
  expect_identical(revcmp("GCCX"), "XGGC")
  # empty input, empty output
  expect_identical(revcmp(""), "")
  # zero-length input, zero-length output
  expect_identical(revcmp(character()), character())
  # NAs are handled
  expect_identical(revcmp(NA), as.character(NA))
})

test_that("revcmp handles IUPAC", {
  expect_identical(revcmp("AYKN"), "NMRT")
})