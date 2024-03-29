% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load_genotypes}
\alias{load_genotypes}
\title{Load table of genotypes}
\usage{
load_genotypes(fp, ...)
}
\arguments{
\item{fp}{path to text file.}

\item{...}{additional arguments passed to [load_csv].}
}
\value{
data frame of genotypes
}
\description{
Load a comma-separated table of genotypes, one pair of alleles per row.  This
information is used to compare samples to genotypes of known individuals in
[summarize_dataset].
}
