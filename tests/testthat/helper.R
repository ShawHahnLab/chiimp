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