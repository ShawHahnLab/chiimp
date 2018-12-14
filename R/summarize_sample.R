#' Versions of summarize_sample
#'
#' This list shows the different versions of \code{\link{summarize_sample}}
#' recognized for use via configuration in full_analysis.  See also
#' \code{\link{sample_analysis_funcs}} and \code{\link{config.defaults}}.
#'
#' @export
sample_summary_funcs <- c("summarize_sample",
                          "summarize_sample_guided")

#' Summarize a processed STR sample
#'
#' Converts an STR sample data frame as produced by
#' \code{\link{analyze_sample}} into a concise list of consistent attributes,
#' suitable for binding together across samples for a dataset.  At this stage
#' the summary is prepared for a single specific locus as in
#' \code{\link{analyze_sample}} but as a list with a fixed length.  The Allele1
#' entries correspond to the sequence with the highest count, Allele2 the
#' second highest.  See the Functions section below for how specific variants
#' of this function behave.
#'
#' @details
#' Entries in the returned list:
#'  * For Allele1 and Allele2:
#'    * \code{Seq}: sequence text for each allele.
#'    * \code{Count}: integer count of occurrences of this exact sequence.
#'    * \code{Length}: integer sequence length.
#'  * \code{Homozygous}: If the sample appears homozygous (if so, the Allele2
#'  entries will be NA).
#'  * \code{Ambiguous}: If a potential allele was ignored due to ambiguous bases
#'  in sequence content (such as "N").
#'  * \code{Stutter}: If a potential allele was ignored due to apparent PCR
#'  stutter.
#'  * \code{Artifact}: If a potential allele was ignored due to apparent PCR
#'  artifact (other than stutter).
#'  * \code{CountTotal}: The total number of sequences in the original sample
#'  data.
#'  * \code{CountLocus}: The number of sequences matching all criteria for the
#'  specified locus in the original sample data.
#'  * \code{ProminentSeqs}: The number of entries above the specified threshold
#'  after all filtering.  This should be either one (for a homozygous sample) or
#'  two (for a heterozygous sample) but conditions such as cross-sample
#'  contamination or excessive PCR stutter can lead to more than two.
#' @md
#'
#' @param sample_data data frame of processed data for one sample as produced by
#'   \code{\link{analyze_sample}}.
#' @param sample.attrs list of sample attributes, such as the rows produced by
#'   \code{\link{prepare_dataset}}.
#' @param counts.min numeric threshold for the minimum number of counts that
#'   must be present, in total across entries passing all filters, for potential
#'   alleles to be considered.
#'
#' @return list of attributes describing the sample.
#'
#' @describeIn summarize_sample Default version of sample summary.
#'
#' @export
summarize_sample <- function(sample_data, sample.attrs, counts.min) {
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(sample_data$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.total <- as.integer(sample_data$Count[1] /
                            sample_data$FractionOfTotal[1])
  if (is.na(count.total))
    count.total <- 0
  count.locus <- sum(sample_data$Count)
  if (count.locus < counts.min) {
    sample_data <- sample_data[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- sample_data[which(sample_data$Category == "Allele")[1], attr.names]
  allele2 <- sample_data[which(sample_data$Category == "Allele")[2], attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  sample.summary <- with(sample_data,
                      c(allele1, allele2,
                        list(Homozygous = sum(Category %in% "Allele") == 1,
                             Ambiguous  = check_category(Category, "Ambiguous"),
                             Stutter    = check_category(Category, "Stutter"),
                             Artifact   = check_category(Category, "Artifact"),
                             CountTotal = count.total,
                             CountLocus = count.locus,
                             ProminentSeqs = prominent.seqs)))
  return(sample.summary)
}

#' @describeIn summarize_sample Summarize a processed STR sample Using known
#'   lengths.  If \code{ExpectedLength1} and optionally \code{ExpectedLength2}
#'   are given in \code{sample.attrs}, the \code{counts.min} threshold is
#'   ignored.  See also \code{\link{analyze_sample_guided}}.
#'
#' @export
summarize_sample_guided <- function(sample_data, sample.attrs, counts.min) {
  expected_lengths <- as.integer(unlist(sample.attrs[c("ExpectedLength1",
                                                       "ExpectedLength2")]))
  expected_lengths <- unique(expected_lengths[! is.na(expected_lengths)])
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(sample_data$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.total <- as.integer(sample_data$Count[1] /
                              sample_data$FractionOfTotal[1])
  if (is.na(count.total))
    count.total <- 0
  count.locus <- sum(sample_data$Count)
  if (length(expected_lengths) == 0 && count.locus < counts.min) {
    sample_data <- sample_data[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- sample_data[which(sample_data$Category == "Allele")[1], attr.names]
  allele2 <- sample_data[which(sample_data$Category == "Allele")[2], attr.names]
  # If the expected lengths were given in the reverse order of the found
  # alleles, flip the found alleles to match.
  if (length(expected_lengths) == 2 &&
      expected_lengths[2] %in% allele1[["Length"]]) {
    allele_tmp <- allele1
    allele1 <- allele2
    allele2 <- allele_tmp
  }
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  sample.summary <- with(sample_data,
                      c(allele1, allele2,
                        list(Homozygous = sum(Category %in% "Allele") == 1,
                             Ambiguous  = check_category(Category, "Ambiguous"),
                             Stutter    = check_category(Category, "Stutter"),
                             Artifact   = check_category(Category, "Artifact"),
                             CountTotal = count.total,
                             CountLocus = count.locus,
                             ProminentSeqs = prominent.seqs)))
  return(sample.summary)
}


# Util --------------------------------------------------------------------


# Did a particular sequence get categorized as a non-allele when it otherwise
# would have been called an allele? Given a sequence category factor and a
# single level, check if that level occurs before a second allele (if any). This
# implies a possible allele was not called due to that category level being
# assigned instead.
check_category <- function(category, lvl) {
  # Stop at the second allele, but if none, check the full vector
  idx_a2 <- which("Allele" == category)[2]
  if (is.na(idx_a2))
    idx_a2 <- length(category)
  # Is there any occurrence of the given lvl in the section to be checked?
  idx <- which(lvl == category[1:idx_a2])[1]
  # Any index found implies TRUE, NA implies FALSE.
  ! is.na(idx)
}

#' Check Sample Data for Potential Allele Matches
#'
#' Check the entries in a processed sample data frame for potential matches to a
#' given locus.
#'
#' @param sample_data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param locus.name character name of locus to match against.
#'
#' @return logical vector of entries for potential alleles.
allele_match <- function(sample_data, locus.name) {
  with(sample_data,
       as.character(MatchingLocus) == locus.name &
         MotifMatch &
         LengthMatch)
}
