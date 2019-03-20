#!/usr/bin/env Rscript

# Lint the package that contains this file's directory, minus some lint
# categories that just annoy me.


# Configure Linters -------------------------------------------------------


# This could be done via a lintr config file if I took the time to figure out
# the syntax.

linters <- list(
  # Linters to add to default list.
  # T_and_F_symbol is quite new as of 2018-12-10.
  yes = c(
    "T_and_F_symbol" # "use TRUE and FALSE, not T and F"
  ),
  # Linters to remove from default list.
  # Some of these are from older versions.
  no = c(
    "multiple_dots", # "Don't use dots in names"
    "camel_case",    # "Don't capitalize stuff"
    "object_name",   # "Don't use dots in names, don't capitalize"
    "object_usage"   # "I don't see that variable"
  )
)


# Detect Package Path -----------------------------------------------------


# If run as a script
args <- commandArgs()
path_script <- gsub("^--file=", "", args[grep("^--file=", args)])
path_script <- normalizePath(path_script)
path_pkg <- dirname(dirname(path_script))

# If run in Rstudio for example
if (! length(path_pkg)) {
  # https://stackoverflow.com/a/16046056
  if (length(sys.frames())) {
    path_pkg <- dirname(dirname(sys.frame(1)$ofile))
  } else {
    # Last fallback, for example run in separate code chunks
    path_pkg <- getwd()
  }
}


# Run lintr ---------------------------------------------------------------

library(lintr)
linters$combined <- lintr::default_linters

# Remove linters
linters$no <- paste0(linters$no, "_linter")
idx <- match(linters$no, names(lintr::default_linters))
idx <- idx[! is.na(idx)]
linters$combined <- lintr::default_linters[-idx]

# Add linters
linters$yes <- paste0(linters$yes, "_linter")
linters$yes <- linters$yes[linters$yes %in% ls("package:lintr")]
.names <- names(linters$yes)
linters$yes <- get(linters$yes, "package:lintr")
names(linters$yes) <- .names
linters$combined <- c(linters$combined, linters$yes)

# Run
results <- lintr::lint_package(path = path_pkg, linters = linters$combined)
results
if (length(path_script) == 1) {
  if (length(results) > 0) {
    quit(status = 1)
  }
}
