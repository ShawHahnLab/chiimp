% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histogram.R
\name{histogram}
\alias{histogram}
\title{Plot histogram of STR sequences}
\usage{
histogram(
  seq_data,
  sample_data = NULL,
  main = "Number of Reads by Sequence Length",
  xlim = range(seq_data$Length),
  cutoff_fraction = NULL,
  render = TRUE
)
}
\arguments{
\item{seq_data}{data frame of dereplicated sequences as created by
\link{analyze_seqs}.}

\item{sample_data}{data frame of filtered and categorized sequences as
created by \link{analyze_sample}.}

\item{main}{title of the plot.}

\item{xlim}{numeric range for x-axis.}

\item{cutoff_fraction}{numeric threshold for the fraction of locus-matching
counts needed to call an allele.  Used to draw a horizontal line if
\code{sample_data} is given.}

\item{render}{Should the plot be drawn to the display device?}
}
\value{
list of data frames for the sets of counts-versus-length bars drawn
in the plot, split by category.
}
\description{
Given processed STR sequences and optionally the filtered version, plot a
histogram of counts per sequence length.
}
