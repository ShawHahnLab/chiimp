# test analyze.dataset ----------------------------------------------------

test_that("analyze.dataset processes samples correctly", {
  # The general case for analyze.dataset
  fail("test not yet implemented")
})

test_that("analyze.dataset runs single-threaded", {
  # analyze.dataset should run properly with num.cores=1 (just slower).
  fail("test not yet implemented")
})

test_that("analyze.dataset warns of missing loci", {
  # If there are locus names in dataset$Locus that are not present in the
  # rownames of locus_attrs, it should throw a warning.
  fail("test not yet implemented")
})
