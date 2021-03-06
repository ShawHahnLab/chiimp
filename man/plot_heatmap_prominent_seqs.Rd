% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{plot_heatmap_prominent_seqs}
\alias{plot_heatmap_prominent_seqs}
\title{Plot heatmap of samples with multiple prominent sequences}
\usage{
plot_heatmap_prominent_seqs(results, ...)
}
\arguments{
\item{results}{combined results list}

\item{...}{additional arguments passed to \code{\link{plot_heatmap}}.}
}
\description{
Given a cross-sample summary data frame as produced by
\code{\link{analyze_dataset}}, plot a heatmap showing samples with more than
two prominent sequences in their analysis output, with sample identifiers on
rows and loci on columns.
}
\seealso{
\code{\link{plot_heatmap}}
}
