% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{save_alignment_images}
\alias{save_alignment_images}
\title{Save alignment visualizations to image files}
\usage{
save_alignment_images(
  alignments,
  dp,
  image.func = "png",
  width = 1600,
  height = 1200,
  res = 150
)
}
\arguments{
\item{alignments}{list of MSA alignment objects, such as created by
\link{summarize_dataset} via \link{align_alleles}.  The name of each alignment will
be used for its filename.}

\item{dp}{output directory path.}

\item{image.func}{name of function to call for saving each image.}

\item{width}{integer width of image.}

\item{height}{integer height of image.}

\item{res}{integer resolution of image in PPI.}
}
\description{
Take a list of alignments, one per locus, and save a plot of each to a
separate image file in a specified directory.  If any of the per-locus
alignment objects is NA it will be skipped.  These are produced by
\link{summarize_dataset} via \link{align_alleles}.
}
