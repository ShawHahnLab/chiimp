% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{make_rownames}
\alias{make_rownames}
\title{Create Row Names for STR Data}
\usage{
make_rownames(data)
}
\arguments{
\item{data}{STR data frame, such as produced by \link{prepare_dataset} or a
summary from \link{summarize_dataset}.}
}
\value{
vector of unique row names
}
\description{
Create unique rownames for the given data frame, using whichever STR-related
metadata columns are available.
}
\seealso{
\link{order_entries}
}
