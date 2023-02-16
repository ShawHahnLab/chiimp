testpth <- function(fname) normalizePath(test_path("data", "configuration", fname))
testrds <- function(fname) readRDS(test_path("data", "configuration", fname))

get_chiimp_opts <- function() {
  chiimp_opts <- options()
  chiimp_opts <- chiimp_opts[grepl("^chiimp\\.", names(chiimp_opts))]
  chiimp_opts
}

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
  # setting a new option (but clean up afterwards)
  cfg("something", 5L)
  expect_identical(getOption("chiimp.something"), 5L)
  cfg("something", NULL)
  expect_false("chiimp.something" %in% names(options()))
})


# apply_config ------------------------------------------------------------


test_that("apply_config loads configuration info into global options", {
  chiimp_opts <- get_chiimp_opts()
  apply_config(CFG_DEFAULTS)
  chiimp_opts_2 <- get_chiimp_opts()
  expect_identical(chiimp_opts, chiimp_opts_2)
  cfg_2 <- CFG_DEFAULTS
  cfg_2$Value[cfg_2$Key == "allele_suffix_len"] <- 12L
  apply_config(cfg_2)
  chiimp_opts_3 <- get_chiimp_opts()
  chiimp_opts_3_exp <- chiimp_opts_2
  chiimp_opts_3_exp$chiimp.allele_suffix_len <- 12L
  expect_identical(chiimp_opts_3, chiimp_opts_3_exp)
})

test_that("apply_config can remove existing global options", {
  cfg_2 <- subset(CFG_DEFAULTS, Key != "allele_suffix_len")
  # by default, it won't remove, so missing entries have no effect on what's
  # already there.
  chiimp_opts <- get_chiimp_opts()
  apply_config(cfg_2)
  chiimp_opts_2 <- get_chiimp_opts()
  expect_identical(chiimp_opts, chiimp_opts_2)
  # With keep=FALSE, missing entries are removed
  apply_config(cfg_2, keep = FALSE)
  chiimp_opts_3 <- get_chiimp_opts()
  chiimp_opts_3_exp <- chiimp_opts_2
  chiimp_opts_3_exp$chiimp.allele_suffix_len <- NULL
  expect_identical(chiimp_opts_3, chiimp_opts_3_exp)
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
  expect_no_warning(config_check_keys(CFG_DEFAULTS))
  alt <- CFG_DEFAULTS
  alt$Key[1] <- "bad"
  expect_warning(
    config_check_keys(alt), "unrecognized config entries")
})


# config_check_version ----------------------------------------------------


test_that("config_check_version warns if config version above pkg version", {
  expect_no_warning(config_check_keys(CFG_DEFAULTS))
  alt <- CFG_DEFAULTS
  alt$Value[match("version", alt$Key)] <- "999.1.2"
  expect_warning(
    config_check_version(alt), "")
})


# parse_config ------------------------------------------------------------


test_that("parse_config loads data structures from a config data frame", {
  config_table <- data.frame(
    Key = c("dataset", "ncores", "report"),
    Value = c("samples.csv", "5", "yes"),
    stringsAsFactors = FALSE)
  cfg <- parse_config(config_table)
  expect_identical(
    cfg, list(dataset = "samples.csv", ncores = 5L, report = TRUE))
})

test_that("parse_config handles uncrecognized keys", {
  config_table <- data.frame(
    Key = c("dataset", "ncores", "do_the_report"),
    Value = c("samples.csv", "5", "yes"),
    stringsAsFactors = FALSE)
  expect_warning(
    cfg <- parse_config(config_table)
  )
  expect_identical(
    cfg, list(dataset = "samples.csv", ncores = 5L, do_the_report = "yes"))
})

test_that("parse_config can use parser specified in input", {
  config_table <- data.frame(
    Key = c("dataset", "ncores", "do_the_report"),
    Value = c("samples.csv", "5", "yes"),
    Parser = c("as.character", "as.integer", "as_bool"),
    stringsAsFactors = FALSE)
  expect_warning(
    cfg <- parse_config(config_table), "unrecognized config entry")
  expect_identical(
    cfg, list(dataset = "samples.csv", ncores = 5L, do_the_report = TRUE))
})

test_that("parse_config handles invalid input values", {
  config_table <- data.frame(
    Key = c("dataset", "ncores"),
    Value = c("samples.csv", "-5"),
    stringsAsFactors = FALSE)
  expect_error(parse_config(config_table), "txt should be an integer")
})


# as_bool -----------------------------------------------------------------


test_that("as_bool parses true/false values from text", {
  trues <- c("TRUE", "T", "true", "yes", "on")
  falses <- c("FALSE", "F", "false", "no", "off")
  for (txt in trues) {
    expect_identical(as_bool(txt), TRUE)
  }
  for (txt in falses) {
    expect_identical(as_bool(txt), FALSE)
  }
  # empty works
  expect_identical(as_bool(""), NA)
  # only length 1
  expect_error(as_locus_vecs(c()),  "txt should be of length 1")
  # controlled vocab
  expect_error(as_bool("NA"))
})


# as_integer_vec ----------------------------------------------------------


test_that("as_integer_vec parses integer vectors from text", {
  expect_equal(as_integer_vec("1+2+3"), 1:3)
  expect_equal(as_integer_vec("1;2@@3"), 1:3)
  expect_equal(as_integer_vec(""), integer())
  expect_equal(as_integer_vec("-5"), -5L)
  expect_equal(as_integer_vec("1.2"), 1:2)
  expect_error(as_locus_vecs(c()),  "txt should be of length 1")
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
  expect_error(as_locus_vecs(c()),  "txt should be of length 1")
})


# as_cpu_cores ------------------------------------------------------------


test_that("as_cpu_cores parses number of CPU cores from text", {
  expect_identical(as_cpu_cores("1"), 1L)
  expect_true(as_cpu_cores("0") > 0)
  expect_true(as_cpu_cores("") > 0)
  expect_error(as_cpu_cores("yes"), "txt should be an integer")
  expect_error(as_cpu_cores("-1"), "txt should be an integer")
  expect_error(as_cpu_cores("5.7"), "txt should be an integer")
  expect_error(as_locus_vecs(c()),  "txt should be of length 1")
})


# is_blank ----------------------------------------------------------------


test_that("is_blank interprets multiple kinds of missing data as TRUE", {
  expect_true(is_blank(NA))
  expect_true(is_blank(NULL))
  expect_true(is_blank(""))
  expect_true(is_blank(c()))
  expect_true(is_blank(c(NA, NA)))
  expect_true(is_blank(c(NA, "")))
  expect_false(is_blank(5))
  expect_false(is_blank(c("", "x")))
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
  apply_config(CFG_DEFAULTS)
  expect_identical(cfg("max_mismatches"), 0L)
})
