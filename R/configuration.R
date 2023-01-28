#' Default Configuration
#'
#' `CFG_DEFAULTS` is CHIIMP's default configuration options and metadata.  The
#' values in this data frame are used unless overridden by custom ones, and the
#' parsing of text values into other types is defined by the `Parser` column
#' here.
#'
#' The columns are:
#'
#' * Key: name of each setting
#' * Value: unparsed text value for each setting
#' * Description: short description of the setting
#' * Example: example value
#' * Parser: name of a function to apply to each entry when setting globally for
#'   the package (see [apply_config])
#'
#' @name CFG_DEFAULTS
#' @export CFG_DEFAULTS
NULL

#' Get/set CHIIMP global configuration options
#'
#' `cfg` is a helper function to provide a shortcut for the equivalent
#' `options(chiimp.key = value)`.  With no arguments, all the currently-defined
#' settings are returned as a list.  With a key given, the corresponding value
#' is returned.  With both key and value, a new option value is set.
#'
#' @param key name of a CHIIMP configuration option
#' @param val value to set for a particular key
#' @export
#' @md
cfg <- function(key = NULL, val = NULL) {
  if (is.null(key)) {
    opts <- options()
    opts <- opts[grepl("^chiimp\\.", names(opts))]
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

#' Apply CHIIMP options globally
#'
#' `apply_config` parses each row in a configuration options data frame and sets
#' a corresponding option with `options(name = value)`.  All options are
#' prefixed with "chiimp." to avoid colliding with options used in other
#' libraries.
#'
#' @param config configuration options data frame, as from [load_config]
#' @param keep keep existing CHIIMP options, or remove them before setting new
#'   ones?
#' @export
#' @md
apply_config <- function(config, keep = TRUE) {
  if (is.data.frame(config)) {
    # work with either list of already-parsed entries, or data frame of text
    # values
    config <- parse_config(config)
  }

  if (! keep) {
    options()
    opts <- options()
    opts <- names(opts)[grepl("^chiimp\\.", names(opts))]
    names(opts) <- opts
    opts <- lapply(opts, function(item) NULL)
    do.call(options, opts)
  }
  names(config) <- paste0("chiimp.", names(config))
  do.call(options, config)
}


# Config Checking ---------------------------------------------------------


#' Compare two version strings
#'
#' `version_compare` takes two strings for software versions(e.g. "1.0.2") and
#' returns `>`, `<`, or `=` depending on whether the first version is later,
#' earlier, or equivalent to the second version.  For example, "1.2.0" >
#' "1.1.5".
#'
#' If there are fewer fields in one string than the other, the one with fewer is
#' padded with zeros for the least-significant positions.  For example 1.2.3 vs
#' 1.2 is equivalent to 1.2.3 vs 1.2.0.
#'
#' @param ver1txt single character string for a software version
#' @param ver2txt single character string for a software version
#' @returns a single character representing if if the first version is later
#'   (`>`), earlier (`<`), or equivalent (`=`) to the second version
#'
#' @md
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


# Config Parsing ----------------------------------------------------------


# take a config table, parse each setting, and return a list of key/value pairs.

#' Parse a configuration data frame into list
#'
#' `parse_config` takes each row of a data frame of configuration options and
#' defines a list item for it, using the Key column for each name and parsing
#' each value according to the Parser column.
#'
#' @param cfg_table data frame of CHIIMP configuration options
#' @returns list with one item per row in the input data frame, with each value
#'   parsed according to the function name in the Parser column
#' @md
parse_config <- function(cfg_table) {
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

check_one_val <- function(txt) {
  if (length(txt) != 1) {
    stop(paste("txt should be of length 1; received", length(txt)))
  }
}

# flexible boolean handling to still act like the old YAML config (well, pre 1.2
# YAML) plus whatever R recognizes
as_bool <- function(txt) {
  check_one_val(txt)
  txt <- toupper(txt)
  map <- c(
    "TRUE" = TRUE, "T" = TRUE, "YES" = TRUE, "ON" = TRUE,
    "FALSE" = FALSE, "F" = FALSE, "NO" = FALSE, "OFF" = FALSE)
  if (! txt %in% c("", names(map))) {
    stop(paste("txt should be TRUE or FALSE; received", txt))
  }
  out <- unname(map[txt])
  out
}

# parse a single character string as a vector of integers separated by
# non-digits
as_integer_vec <- function(txt) {
  check_one_val(txt)
  as.integer(strsplit(txt, "[^-0-9]+")[[1]])
}

# parse a single character string as a named list of vectors
# e.g. "name=item1/item2/item3;name2=item4/item5;..."
as_locus_vecs <- function(txt) {
  check_one_val(txt)
  chunks <- strsplit(txt, "; *")[[1]]
  chunk_names <- sub("=.*", "", chunks)
  vecs <- lapply(chunks, function(chunk) {
    strsplit(sub(".*=", "", chunk), "/")[[1]]
  })
  names(vecs) <- chunk_names
  vecs
}

# autodetect CPU cores if 0
as_cpu_cores <- function(txt) {
  check_one_val(txt)
  val <- as.integer(txt)
  if (val == 0) {
    val <- max(1, as.integer(parallel::detectCores() / 2) - 1)
  }
  val
}

# make designated paths absolute relative to working directory, unless already
# absolute
as_abs_path <- function(txt) {
  check_one_val(txt)
  if (substr(txt, 1, 1) != .Platform$file.sep) {
    file.path(normalizePath("."), txt)
  } else {
    txt
  }
}

as_rel_path <- function(txt) {
  check_one_val(txt)
  txt
}