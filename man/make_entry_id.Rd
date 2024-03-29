% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{make_entry_id}
\alias{make_entry_id}
\title{Create identifiers for STR Data}
\usage{
make_entry_id(data)
}
\arguments{
\item{data}{STR data frame, such as produced by \link{prepare_dataset} or a
summary from \link{summarize_dataset}.}
}
\value{
character vector of entry identifiers
}
\description{
Create entry IDs for the given data frame, using whichever STR-related
metadata columns are available.  These will not necessarily be unique.
}
