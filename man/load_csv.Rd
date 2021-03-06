% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load_csv}
\alias{load_csv}
\alias{save_csv}
\title{Load and save tables from CSV}
\usage{
load_csv(fp, ...)

save_csv(data, fp, ...)
}
\arguments{
\item{fp}{path to text file.}

\item{...}{additional arguments passed to \code{\link[utils]{read.table}} or
\code{\link[utils]{write.table}}.}

\item{data}{data frame to save to CSV file}
}
\value{
data frame
}
\description{
Load/save a comma-separated table from/to a CSV file.  (These are generic
wrapper functions used by more specific loaders like
\code{\link{load_locus_attrs}}.)
}
\section{Functions}{
\itemize{
\item \code{load_csv}: Load CSV

\item \code{save_csv}: Save CSV
}}

