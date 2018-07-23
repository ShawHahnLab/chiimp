#!/usr/bin/env Rscript

# Lint the package that contains this file's directory, minus some lint
# categories that just annoy me.

args <- commandArgs()
f <- gsub("^--file=", "", args[grep("^--file=", args)])
f <- normalizePath(f)
path <- dirname(dirname(f))

linters_no <- c("multiple_dots", # "Don't use dots in names"
                "camel_case",    # "Don't capitalize stuff"
                "object_usage")  # "I don't see that variable"
linters_no <- paste0(linters_no, "_linter")
linters <- lintr::default_linters[-match(linters_no,
                                         names(lintr::default_linters))]
lintr::lint_package(path = path, linters = linters)
