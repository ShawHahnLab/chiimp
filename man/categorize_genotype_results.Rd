% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorize.R
\name{categorize_genotype_results}
\alias{categorize_genotype_results}
\title{Categorize genotyping results}
\usage{
categorize_genotype_results(results_summary)
}
\arguments{
\item{results_summary}{cross-sample summary data frame as produced by
\link{analyze_dataset} with extra columns as produced by
\link{match_known_genotypes}.}
}
\value{
factor defining genotyping result category for every row of the input
data frame.
}
\description{
For a given results summary data frame that has \code{CorrectAllele1Seq} and
\code{CorrectAllele2Seq} columns (such as produced by \link{match_known_genotypes})
added, create a factor labeling every row of the input data frame by its
genotyping outcome.
}
\details{
Levels in the returned factor, in order:
\itemize{
\item Correct: one/two alleles match.
\item Incorrect at least one allele does not match.
\item Blank: No alleles were called in the analysis even though known genotypes
were supplied.
\item Dropped Allele: One called allele is correct for a heterozygous individual,
but no second allele was called.
}

Cases that should not occur, such as \code{CorrectAllele1Seq} and
\code{CorrectAllele2Seq} both set to NA, map to NA in the returned factor.
}
