testpth <- function(fname) normalizePath(test_path("data", "configuration", fname))
testrds <- function(fname) readRDS(test_path("data", "configuration", fname))


# cfg ---------------------------------------------------------------------


test_that("cfg handles chiimp options", {
  # check one basic case
  val <- as.integer(CFG_DEFAULTS$Value[
    match("max_mismatches", CFG_DEFAULTS$Key)])
  expect_equal(cfg("max_mismatches"), val)
  expect_equal(getOption("chiimp.max_mismatches"), val)
  # we can go full circle with all options
  chiimp_opts <- cfg()
  expect_true(is.list(chiimp_opts))
  for (item in names(chiimp_opts)) {
    cfg(item, chiimp_opts[[item]])
  }
  chiimp_opts_2 <- cfg()
  expect_identical(chiimp_opts, chiimp_opts_2)
  # setting a new option
  cfg("something", 5L)
  expect_identical(getOption("chiimp.something"), 5L)
})


# apply_config ------------------------------------------------------------


test_that("apply_config loads configuration info into global options", {

})

test_that("apply_config can remove existing global options", {

})


# version_compare ---------------------------------------------------------


test_that("version_compare can order two software version strings", {
  expect_equal(version_compare("1.2.3", "1.4.0"), "<")
  expect_equal(version_compare("1.4.0", "1.2.3"), ">")
  expect_equal(version_compare("1.4.0", "1.4.0"), "=")
  expect_equal(version_compare("0.5.0", "0.5.1"), "<")
  expect_equal(version_compare("1.2.3", "1"), ">") # treats 1 like 1.0.0
  expect_equal(version_compare("", ""), "=")
})


# config_check_keys -------------------------------------------------------


test_that("config_check_keys warns of unrecognized config keys", {

})


# config_check_version ----------------------------------------------------


test_that("config_check_version warns if config version above pkg version", {

})


# parse_config ------------------------------------------------------------


test_that("parse_config loads a config data frame from CSV", {

})



# as_bool -----------------------------------------------------------------


test_that("as_bool parses true/false values from text", {

})


# as_integer_vec ----------------------------------------------------------


test_that("as_integer_vec parses integer vectors from text", {
  expect_equal(as_integer_vec("1+2+3"), 1:3)
  expect_equal(as_integer_vec("1;2@@3"), 1:3)
  expect_equal(as_integer_vec(""), integer())
  expect_equal(as_integer_vec("-5"), -5L)
  expect_equal(as_integer_vec("1.2"), 1:2)
  # only length 1
  expect_error(as_integer_vec(c("1", "2")))
  expect_error(as_integer_vec(c()))
})


# as_locus_vecs -----------------------------------------------------------


test_that("as_locus_vecs parses lists of locus names from text", {
  expect_equal(
    as_locus_vecs("X=1/2/3;Y=A/B"),
    list(X = c("1", "2", "3"), Y = c("A", "B")))
  expect_equal(
    as_locus_vecs("X=1/2/3; Y=A/B"),
    list(X = c("1", "2", "3"), Y = c("A", "B")))
  output <- list()
  # evidently "list()" is a different thing than "named list()"
  names(output) <- character()
  expect_equal(as_locus_vecs(""), output)
  # only length 1
  expect_error(as_locus_vecs(c("X=1/2/3", "Y=A/B")))
  expect_error(as_locus_vecs(c()))
})


# as_cpu_cores ------------------------------------------------------------


test_that("as_cpu_cores parses number of CPU cores from text", {

})


# as_abs_path -------------------------------------------------------------


test_that("as_abs_path parses abolute path from text", {

})


# as_rel_path -------------------------------------------------------------


test_that("as_rel_path parses abolute path from text", {

})


# is_blank ----------------------------------------------------------------


test_that("is_blank interprets multiple kinds of missing data as TRUE", {

})


# package configuration ---------------------------------------------------


test_that("CHIIMP package config updates work", {
  # we should be able to update the config "live" and have it take effect for
  # the whole package.  using find_primer_matches as an arbitrary example case.
  expect_identical(cfg("max_mismatches"), 0L)
  result <- find_primer_matches("TAAGAAATGCTTATATGGCCATAAATCAAC", "TAAGAAA")
  expect_identical(result$Mismatches, 0L)
  result <- find_primer_matches("TAAGAAATGCTTATATGGCCATAAATCAAC", "TAAGAAT")
  expect_identical(result$Mismatches, integer())
  cfg("max_mismatches", 1)
  result <- find_primer_matches("TAAGAAATGCTTATATGGCCATAAATCAAC", "TAAGAAA")
  expect_identical(result$Mismatches, 0L)
  result <- find_primer_matches("TAAGAAATGCTTATATGGCCATAAATCAAC", "TAAGAAT")
  expect_identical(result$Mismatches, 1L)
  # can reset to get the defaults again
  setup_package()
  expect_identical(cfg("max_mismatches"), 0L)
})
