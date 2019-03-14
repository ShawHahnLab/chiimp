# Install CHIIMP on Windows.

# Find the path to the directory containing this script.  We need this for
# package testing and installation below.
args <- commandArgs()
f <- gsub("^--file=", "", args[grep("^--file=", args)])
f <- normalizePath(f)
path <- dirname(f)

UPROF <- Sys.getenv("USERPROFILE")

# If no library paths are writeable, try creating a user library.
if (! any(file.access(.libPaths(), 2) == 0)) {
  # This is the directory I see RStudio create automatically on first start,
  # and the command-line R also detects it.
  ver <- paste(version$major, sub("\\..*", "", version$minor), sep = ".")
  dp <- file.path(UPROF, "Documents", "R", "win-library", ver)
  dir.create(dp, recursive = TRUE)
  # On a second run through this will get picked up automatically,
  # but if we want it right now we have to add it to the list manually.
  .libPaths(dp)
}

haspkg <- function(pkgname) {
  suppressMessages(suppressWarnings(
    require(pkgname, character.only = TRUE, quietly = TRUE)
  ))
}

if (! haspkg("devtools")) {
  cat("\n")
  cat("### Installing devtools\n")
  cat("\n")
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

if (! haspkg("msa")) {
  cat("\n")
  cat("### Installing Bioconductor and MSA\n")
  cat("\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("msa", suppressUpdates = TRUE)
}

cat("\n")
cat("### Installing CHIIMP\n")
cat("\n")
devtools::install(path, upgrade = "never")

shortcut_path <- file.path(UPROF, "Desktop", "CHIIMP.lnk")
chiimp_path <- system.file("bin", "chiimp.cmd", package = "chiimp")
# https://stackoverflow.com/a/30029955/6073858
args <- c(paste0("$s=(New-Object -COM WScript.Shell).CreateShortcut('",
				         shortcut_path,
				         "');"),
		      paste0("$s.TargetPath='", chiimp_path, "';"),
	        "$s.Save();")
system2("powershell", args)
