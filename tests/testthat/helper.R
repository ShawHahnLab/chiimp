# execute expression inside of a temporary directory, and remove the directory
# afterward
within_tmpdir <- function(expr) {
  here <- getwd()
  data.dir <- tempfile()
  dir.create(data.dir)
  setwd(data.dir)
  tryCatch(eval(expr), finally = {
    unlink(x = data.dir, recursive = TRUE)
    setwd(here)
  })
}

# append an empty string to each the given files
touch <- function(fps) {
  lapply(fps, function(fp) cat("", file = fp, append = TRUE))
}

# seq_sets is per-sample list of per-locus sequences
write_seqs <- function(seq_sets, outdir, fmt="%s-%s.fasta") {
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