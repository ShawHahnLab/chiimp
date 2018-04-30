# Create a summary list for a single sample, targeting a single locus.  See
# summarize_sample below for the entries in the returned list.

# Different versions of summarize_sample to allow via configuration in
# full_analysis.
sample_summary_funcs <- c("summarize_sample",
                          "summarize_sample_guided",
                          "summarize_sample_naive")

#' Summarize a processed STR sample
#'
#' Converts a full STR sample data frame into a concise list of consistent
#' attributes, suitable for binding together across samples for a dataset.  At
#' this stage the summary is prepared for a single specific locus, in contrast
#' to \code{\link{analyze_seqs}}.  The Allele1 entries correspond to the
#' sequence with the highest count, Allele2 the second highest.
#'
#' @details
#' Entries in the returned list:
#'  * For Allele1 and Allele2:
#'    * Seq: sequence text for each allele.
#'    * Count: integer count of occrrences of this exact sequence.
#'    * Length: integer sequence length.
#'  * Homozygous: If the sample appears homozygous (if so, the Allele2 entries
#'  will be NA).
#'  * Ambiguous: If a potential allele was ignored due to ambiguous bases in
#'  sequence content (such as "N").
#'  * Stutter: If a potential allele was ignored due to apparent PCR stutter.
#'  * Artifact: If a potential allele was ignored due to apparent PCR artifact
#'  (other than stutter).
#'  * CountTotal: The total number of sequences in the original sample data.
#'  * CountLocus: The number of sequences matching all criteria for the
#'  specified locus in the original sample data.
#'  * ProminentSeqs: The number of entries above the specified threshold after
#'  all filtering.  This should be either one (for a homozygous sample) or two
#'  (for a heterozygous sample) but conditions such as cross-sample
#'  contamination or excessive PCR stutter can lead to more than two.
#' @md
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param sample.attrs list of sample attributes, such as the rows produced by
#'   \code{\link{prepare_dataset}}.  Used to select the locus name to filter on.
#' @param fraction.min numeric threshold for the minimum proportion of counts a
#'   given entry must have, compared to the total matching all criteria for that
#'   locus, to be considered as a potential allele.
#' @param counts.min numeric threshold for the minimum number of counts that
#'   must be present, in total across entries passing all filters, for potential
#'   alleles to be considered.
#'
#' @return list of attributes describing the sample.
#'
#' @export
summarize_sample <- function(sample.data, sample.attrs, fraction.min,
                             counts.min) {
  chunk <- analyze_sample(sample.data, sample.attrs, fraction.min)
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(chunk$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.locus <- sum(chunk$Count)
  if (count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[which(chunk$Category == "Allele")[1], attr.names]
  allele2 <- chunk[which(chunk$Category == "Allele")[2], attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  sample.summary <- with(chunk,
                      c(allele1, allele2,
                        list(Homozygous = sum(Category %in% "Allele") == 1,
                             Ambiguous  = check_category(Category, "Ambiguous"),
                             Stutter    = check_category(Category, "Stutter"),
                             Artifact   = check_category(Category, "Artifact"),
                             CountTotal = sum(sample.data$Count),
                             CountLocus = count.locus,
                             ProminentSeqs = prominent.seqs)))
  return(sample.summary)
}

#' Summarize a processed STR sample Using Known Lengths
#'
#' Converts a full STR sample data frame into a concise list of consistent
#' attributes, suitable for binding together across samples for a dataset.  At
#' this stage the summary is prepared for a single specific locus, in contrast
#' to \code{\link{analyze_seqs}}.  The Allele1 entries correspond to the
#' sequence with the highest count, Allele2 the second highest.
#'
#' @details
#' Entries in the returned list:
#'  * For Allele1 and Allele2:
#'    * Seq: sequence text for each allele.
#'    * Count: integer count of occrrences of this exact sequence.
#'    * Length: integer sequence length.
#'  * Homozygous: If the sample appears homozygous (if so, the Allele2 entries
#'  will be NA).
#'  * Ambiguous: If a potential allele was ignored due to ambiguous bases in
#'  sequence content (such as "N").
#'  * Stutter: If a potential allele was ignored due to apparent PCR stutter.
#'  * Artifact: If a potential allele was ignored due to apparent PCR artifact
#'  (other than stutter).
#'  * CountTotal: The total number of sequences in the original sample data.
#'  * CountLocus: The number of sequences matching all criteria for the
#'  specified locus in the original sample data.
#'  * ProminentSeqs: The number of entries above the specified threshold after
#'  all filtering.  This should be either one (for a homozygous sample) or two
#'  (for a heterozygous sample) but conditions such as cross-sample
#'  contamination or excessive PCR stutter can lead to more than two.
#' @md
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param sample.attrs list of sample attributes, such as the rows produced by
#'   \code{\link{prepare_dataset}}.  Used to select the locus name to filter on
#'   and the sequence lengths to select.  If the ExpectedLength1 and
#'   ExpectedLength2 columns are given, \code{counts.min} and the stutter
#'   filtering are ignored. \code{fraction.min} is also ignored if two lengths
#'   are given.
#' @param fraction.min numeric threshold for the minimum proportion of counts a
#'   given entry must have, compared to the total matching all criteria for that
#'   locus, to be considered as a potential allele.
#' @param counts.min numeric threshold for the minimum number of counts that
#'   must be present, in total across entries passing all filters, for potential
#'   alleles to be considered.
#'
#' @return list of attributes describing the sample.
#'
#' @export
summarize_sample_guided <- function(sample.data, sample.attrs, fraction.min,
                             counts.min) {
  expected_lengths <- as.integer(unlist(sample.attrs[c("ExpectedLength1",
                                                       "ExpectedLength2")]))
  expected_lengths <- unique(expected_lengths[! is.na(expected_lengths)])

  chunk <- analyze_sample_guided(sample.data, sample.attrs, fraction.min)

  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(chunk$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.locus <- sum(chunk$Count)
  if (length(expected_lengths) == 0 && count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[which(chunk$Category == "Allele")[1], attr.names]
  allele2 <- chunk[which(chunk$Category == "Allele")[2], attr.names]
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
  sample.summary <- with(chunk,
                      c(allele1, allele2,
                        list(Homozygous = sum(Category %in% "Allele") == 1,
                             Ambiguous  = check_category(Category, "Ambiguous"),
                             Stutter    = check_category(Category, "Stutter"),
                             Artifact   = check_category(Category, "Artifact"),
                             CountTotal = sum(sample.data$Count),
                             CountLocus = count.locus,
                             ProminentSeqs = prominent.seqs)))
  return(sample.summary)
}

#' Summarize a processed STR sample, Simple Version
#'
#' Summarize as in \code{\link{summarize_sample}}, but skip allele_match and
#' stutter/artifact removal.  Entries in the returned list are the same.
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param sample.attrs list of sample attributes, such as the rows produced by
#'   \code{\link{prepare_dataset}}.  Used to select the locus name to filter on.
#' @param fraction.min numeric threshold for the minimum proportion of counts a
#'   given entry must have, compared to the total matching all criteria for that
#'   locus, to be considered as a potential allele.
#' @param counts.min numeric threshold for the minimum number of counts that
#'   must be present, in total across entries passing all filters, for potential
#'   alleles to be considered.
#'
#' @return list of attributes describing the sample.
summarize_sample_naive <- function(sample.data, sample.attrs, fraction.min,
                                   counts.min) {
  chunk <- analyze_sample_naive(sample.data, sample.attrs, fraction.min)
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(chunk$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.locus <- sum(chunk$Count)
  if (count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[which(chunk$Category == "Allele")[1], attr.names]
  allele2 <- chunk[which(chunk$Category == "Allele")[2], attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  sample.summary <- with(chunk,
                      c(allele1, allele2,
                        list(Homozygous = sum(Category %in% "Allele") == 1,
                             Ambiguous  = check_category(Category, "Ambiguous"),
                             Stutter    = check_category(Category, "Stutter"),
                             Artifact   = check_category(Category, "Artifact"),
                             CountTotal = sum(sample.data$Count),
                             CountLocus = count.locus,
                             ProminentSeqs = prominent.seqs)))
  return(sample.summary)
}

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
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param locus.name character name of locus to match against.
#'
#' @return logical vector of entries for potential alleles.
allele_match <- function(sample.data, locus.name) {
  with(sample.data,
       as.character(MatchingLocus) == locus.name &
         MotifMatch &
         LengthMatch)
}
