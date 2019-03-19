# Install CHIIMP package and desktop icon.

# Should add a log file but R makes it tricky to safely capture stderr.  Maybe
# just grab messages/warnings/errors specifically.
# https://stackoverflow.com/questions/45036224/how-to-write-errors-and-warnings-to-a-log-file
# https://stackoverflow.com/questions/19433848/handling-errors-before-warnings-in-trycatch

# Functions ---------------------------------------------------------------


# Find the path to the directory containing this script.
get_script_path <- function() {
  args <- commandArgs()
  f <- gsub("^--file=", "", args[grep("^--file=", args)])
  f <- normalizePath(f)
  d <- dirname(f)
  if (length(d)) {
    return(d)
  } else {
    return(NULL)
  }
}

# Check very quietly if a package is installed
haspkg <- function(pkgname) {
  suppressMessages(suppressWarnings({
    status <- requireNamespace(pkgname, quietly = TRUE)
    return(status) # (in two steps for for visibility)
  }))
}

# Get operating system: "linux", "osx", "windows", ...
# Based off of:
# https://www.r-bloggers.com/identifying-the-os-from-r/
# (But I see no mention of Sys.info() not being implemented in the docs as of R
# 3.2.3 and 3.5.1, and we only expect one of three possibilities, so we can
# simplify a bit.)
get_os <- function() {
  os <- Sys.info()[["sysname"]]
  if (os == "Darwin") {
    os <- "osx"
  }
  tolower(os)
}


# Setup R User Library ----------------------------------------------------


# Sys.getenv("R_LIBS_USER") seems to work even if that variable was not set
# going into R and even if the directory doesn't yet exist.
#
# Documentation is a bit vague but hints in this direction:
#
# ?.libPaths
# "By default R_LIBS is unset, and R_LIBS_USER is set to directory
# ‘R/R.version$platform-library/x.y’ of the home directory (or
# ‘Library/R/x.y/library’ for CRAN macOS builds), for R x.y.z."
#
# The behavior I actually see is:
#  Linux: ~/R/R.version$platform-library/x.y
#  OSX: ~/Library/R/x.y/library
#  Windows: %USERPROFILE%\Documents\R\win-library\x.y

# Check if we have write access to any library paths.  If not, create the user
# library appropriate for the detected operating system.
setup_user_library <- function() {
  if (! any(file.access(.libPaths(), 2) == 0)) {
    dp <- normalizePath(Sys.getenv("R_LIBS_USER"), mustWork = FALSE)
    dir.create(dp, recursive = TRUE)
    # On a second run through this will get picked up automatically,
    # but if we want it right now we have to add it to the list manually.
    .libPaths(dp)
    return(dp)
  }
}


# Setup Desktop Icon ------------------------------------------------------


# Create a Desktop icon appropriate for the detected operating system.
#   Windows: A .lnk Shortcut to a .cmd wrapper script
#   OS X: A symbolic link to a small .app wrapper application
#   Linux: A .desktop INI file pointing to a .sh wrapper script
setup_icon <- function() {
  os <- get_os()
  if (os == "windows") {
    setup_icon_windows()
  } else if (os == "osx") {
    setup_icon_osx()
  } else if (os == "linux") {
    setup_icon_linux()
  } else {
    warning("Operating system not recognized; skipping Desktop icon setup")
  }
}

setup_icon_linux <- function() {
  chiimp_path <- system.file("bin", "chiimp.sh", package = "chiimp")
  desktop_path <- normalizePath("~/Desktop", mustWork = FALSE)
  icon_path <- NULL
  if (dir.exists(desktop_path)) {
    desktop_file <- paste(
      "[Desktop Entry]",
      "Type=Application",
      "Terminal=true",
      "Name=CHIIMP",
      paste("Exec", chiimp_path, sep = ":"),
      sep = "\n")
    icon_path <- file.path(desktop_path, "CHIIMP.desktop")
    icon_path <- normalizePath(icon_path, mustWork = FALSE)
    cat(desktop_file, file = icon_path)
    # TODO double-check if .desktop files actually need to be marked
    # exectuable.  This may not be necessary.
    system2("chmod", args = c("+x", icon_path))
  }
  return(icon_path)
}

setup_icon_osx <- function() {
  chiimp_path <- system.file("bin", "chiimp.app", package = "chiimp")
  droplet_path <- file.path(chiimp_path, "Contents", "MacOS", "droplet.gz")
  desktop_path <- normalizePath("~/Desktop", mustWork = FALSE)
  system2("gunzip", args = droplet_path)
  icon_path <- NULL
  if (dir.exists(desktop_path)) {
    icon_path <- file.path(desktop_path, "CHIIMP")
    icon_path <- normalizePath(icon_path, mustWork = FALSE)
    # Create the symbolic link, replacing whatever link may have already been
    # there.
    # Also consider R's file.symlink().
    system2("ln", args = c("-shf", chiimp_path, icon_path))
  }
  return(icon_path)
}

setup_icon_windows <- function() {
  chiimp_path <- system.file("bin", "chiimp.cmd", package = "chiimp")
  uprof <- Sys.getenv("USERPROFILE")
  desktop_path <- normalizePath(file.path(uprof, "Desktop"), mustWork = FALSE)
  icon_path <- NULL
  if (dir.exists(desktop_path)) {
    icon_path <- file.path(uprof, "Desktop", "CHIIMP.lnk")
    icon_path <- normalizePath(icon_path, mustWork = FALSE)
    # https://stackoverflow.com/a/30029955/6073858
    args <- c(paste0("$s=(New-Object -COM WScript.Shell).CreateShortcut('",
                     icon_path,
                     "');"),
              paste0("$s.TargetPath='", chiimp_path, "';"),
              "$s.Save();")
    system2("powershell", args = args)
  }
  return(icon_path)
}


# Install -----------------------------------------------------------------


install <- function(path_package = NULL) {

  results <- list()

  if (is.null(path_package)) {
    path_installer <- get_script_path()
    if (is.null(path_installer)) {
      stop("install must be run as script or given explicit package path")
    }
    path_package <- normalizePath(file.path(path_installer, "..", ".."))
  }

  results$new_user_library <- setup_user_library()

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
  devtools::install(path_package, upgrade = "never")
  results$icon <- setup_icon()

  invisible(results)
}

# e.g., if __name__ == "__main__"
if (! is.null(get_script_path())) {
 install()
}
