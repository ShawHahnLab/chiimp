# execute expression inside of a temporary directory, and remove the directory
# afterward
within_tmpdir <- function(expr) {
  here <- getwd()
  data_dir <- tempfile()
  dir.create(data_dir)
  setwd(data_dir)
  tryCatch(eval(expr), finally = {
    unlink(x = data_dir, recursive = TRUE)
    setwd(here)
  })
}

# append an empty string to each the given files
touch <- function(fps) {
  lapply(fps, function(fp) cat("", file = fp, append = TRUE))
}

# seq_sets is per-sample list of per-locus sequences
write_seqs <- function(seq_sets, outdir, fmt = "%s-%s.fasta") {
  if (! dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)
  for (sn in names(seq_sets)) {
    for (ln in names(seq_sets[[sn]])) {
      fp <- file.path(outdir, sprintf(fmt, sn, ln))
      n <- names(seq_sets[[sn]][[ln]])
      if (is.null(n))
        n <- seq_along(seq_sets[[sn]][[ln]])
      dnar::write.fa(names = n,
                     dna = seq_sets[[sn]][[ln]],
                     fileName = fp)
    }
  }
}

# for testthat pre-3.1.5
if (! exists("expect_no_warning", envir = as.environment("package:testthat"))) {
  expect_no_warning <- function(...) {
    wrns <- testthat::capture_warnings(...)
    testthat::expect_identical(wrns, character())
  }
}
