% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_dataset.R
\name{find_closest_matches}
\alias{find_closest_matches}
\title{Find closest matches in distance matrix}
\usage{
find_closest_matches(dist_mat, range = 2, maximum = 8)
}
\arguments{
\item{dist_mat}{matrix of distance values, such as produced by
\link{make_dist_mat} and \link{make_dist_mat_known}.}

\item{range}{optional numeric for distances to each set of nearby names,
relative to the closest match.}

\item{maximum}{optional numeric maximum value for any distance.}
}
\value{
list of named vectors containing distances for each sample.
}
\description{
Given a distance matrix with samples on rows and names on columns, return a
list of vectors with the closest-matching names for each sample.
}
