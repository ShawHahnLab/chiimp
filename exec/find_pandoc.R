#!/usr/bin/env Rscript

# Figure out the path to the pandoc executable.  In order of priority:
#
#  1. Anything set in RSTUDIO_PANDOC env variable
#  2. Anything found in installed copies of RStudio
#  3. Any other pandoc found on PATH

find_pandoc <- function() {
  # If there's already an environment variable defined, we'll just give that.
  pandoc_env <- Sys.getenv("RSTUDIO_PANDOC")
  if (pandoc_env != "") {
    return(pandoc_env)
  }

  # On Windows this should give something like "c:/", but an empty string
  # everywhere else.
  prefix <- sub("/.*$", "/", normalizePath(".", winslash = "/"))

  # Search for RStudio installations
  search_paths <- c(
    list.files(list.files(prefix, pattern = "Program Files.*", full.names = TRUE), pattern = "RStudio.*", full.names = TRUE),
    list.files("/Applications", pattern = "^RStudio.*.app", full.names = TRUE),
    "/usr/lib/rstudio")

  execs <- file.path(
    strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]],
    c("pandoc", "pandoc.exe"))
  execs <- execs[file.exists(execs)]

  # Find files (no dirs) that are named exactly "pandoc", and get the full path
  # to the parent dir of any found
  pandoc_dirs_available <- c(
    dirname(list.files(search_paths, recursive = TRUE, full.names = TRUE, pattern = "^pandoc(?:\\.exe)$")),
    dirname(execs))

  # Let rmarkdown package decide which pandoc to use, if that feature is
  # available (>=2.2).  It will also implicitly accept a pandoc executable found
  # on the PATH.  If that function's not available we'll just pick a pandoc.
  # https://stackoverflow.com/a/24692264
  pandoc_dir <- if ("find_pandoc" %in% getNamespaceExports("rmarkdown")) {
    rmarkdown::find_pandoc(dir = pandoc_dirs_available)$dir
  } else {
    file.path(pandoc_dirs_available[1], "pandoc")
  }
  return(pandoc_dir)
}

cat(find_pandoc(), end = "\n")
