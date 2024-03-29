% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{plot_heatmap_stutter}
\alias{plot_heatmap_stutter}
\title{Plot heatmap of suspected PCR stutter}
\usage{
plot_heatmap_stutter(results, ...)
}
\arguments{
\item{results}{combined results list}

\item{...}{additional arguments passed to \link{plot_heatmap}.}
}
\description{
Given a cross-sample summary data frame as produced by \link{analyze_dataset},
plot a heatmap showing which samples had alleles ignored due to suspected PCR
stutter, with sample identifiers on rows and loci on columns.
}
\seealso{
\link{plot_heatmap}
}
