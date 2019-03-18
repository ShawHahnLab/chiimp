# Install CHIIMP package and desktop icon.


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


# Check if we have write access to any library paths.  If not, create the user
# library.

setup_user_library <- function() {
  os <- get_os()
  if (os == "windows") {
    setup_user_library_windows()
  } else if (os == "osx") {
    setup_user_library_osx()
  } else if (os == "linux") {
    setup_user_library_linux()
  } else {
    warning("Operating system not recognized; skipping user library setup")
  }
}

setup_user_library_linux <- function() {
  # TODO implement as for Windows
}

setup_user_library_osx <- function() {
  # TODO implement as for Windows
}

setup_user_library_windows <- function() {
  if (! any(file.access(.libPaths(), 2) == 0)) {
    dp <- get_user_library_windows()
    dir.create(dp, recursive = TRUE)
    # On a second run through this will get picked up automatically,
    # but if we want it right now we have to add it to the list manually.
    .libPaths(dp)
    return(dp)
  }
}

get_user_library_windows <- function() {
  # This is the directory I see RStudio create automatically on first start,
  # and the command-line R also detects it.
  uprof <- Sys.getenv("USERPROFILE")
  ver <- paste(version$major, sub("\\..*", "", version$minor), sep = ".")
  dp <- file.path(uprof, "Documents", "R", "win-library", ver)
  dp <- normalizePath(dp, mustWork = FALSE)
  return(dp)
}


# Setup Desktop Icon ------------------------------------------------------


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
    system2("chmod", "+x", icon_path)
  }
  return(icon_path)
}

setup_icon_osx <- function() {
  chiimp_path <- system.file("bin", "chiimp.app", package = "chiimp")
  droplet_path <- file.path(chiimp_path, "Contents", "MacOS", "droplet.gz")
  desktop_path <- normalizePath("~/Desktop", mustWork = FALSE)
  system2("gunzip", droplet_path)
  icon_path <- NULL
  if (dir.exists(desktop_path)) {
    icon_path <- file.path(desktop_path, "CHIIMP")
    icon_path <- normalizePath(icon_path, mustWork = FALSE)
    system2("ln", "-s", chiimp_path, icon_path)
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
    system2("powershell", args)
  }
  return(icon_path)
}


# Install -----------------------------------------------------------------


install <- function(path_package = NULL) {

  if (is.null(path_package)) {
    path_installer <- get_script_path()
    if (is.null(path_installer)) {
      stop("install must be run as script or given explicit package path")
    }
    path_package <- normalizePath(file.path(path_installer, "..", ".."))
  }

  setup_user_library()

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
  setup_icon()
}

# e.g., if __name__ == "__main__"
if (! is.null(get_script_path())) {
 install()
}
