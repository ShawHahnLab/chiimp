# Create a summary list for a single sample, targeting a single locus.  See
# summarize_sample below for the entries in the returned list.

# Different versions of summarize_sample to allow via configuration in
# full_analysis.
sample_summary_funcs <- c("summarize_sample",
                          "summarize_sample_guided",
                          "summarize_sample_naive",
                          "summarize_sample_by_length")

#' Summarize a processed STR sample
#'
#' Converts a full STR sample data frame into a concise list of consistent
#' attributes, suitable for binding together across samples for a dataset.  At
#' this stage the summary is prepared for a single specific locus, in contrast
#' to \code{\link{analyze_sample}}.  The Allele1 entries correspond to the
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
#'   \code{\link{analyze_sample}}.
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
  # extract sample data entries that meet all criteria for a potential allele
  locus.name <- sample.attrs[["Locus"]]
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  # Note that counts.locus is more restrictive than the total counts of all
  # entries with MatchingLocus equal to the given locus name, since idx includes
  # the extra checks above.
  count.total <- sum(sample.data$Count)
  count.locus <- sum(chunk$Count)
  # Threshold potential alleles at minimum count
  chunk <- chunk[chunk$Count >= fraction.min * count.locus, ]
  # Remove ambiguous sequences, if present.
  ambig <- chunk[2, "Ambiguous"]
  chunk <- chunk[! chunk[, "Ambiguous"], ]
  # Remove stutter, if present.
  stutter <- !is.na(chunk[2, "Stutter"])
  chunk <- chunk[is.na(chunk[, "Stutter"]), ]
  # Remove artifact-like sequences, if present.
  artifact <- !is.na(chunk[2, "Artifact"])
  chunk <- chunk[is.na(chunk[, "Artifact"]), ]

  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- nrow(chunk)
  # enforce count limit after all filtering (but before stutter removal or
  # fraction thresholding)
  if (count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[1, attr.names]
  allele2 <- chunk[2, attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  homozygous <- nrow(chunk) == 1
  sample.summary <- c(allele1, allele2,
                      list(Homozygous = homozygous,
                           Ambiguous = ambig,
                           Stutter = stutter,
                           Artifact = artifact,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

#' Summarize a processed STR sample Using Known Lengths
#'
#' Converts a full STR sample data frame into a concise list of consistent
#' attributes, suitable for binding together across samples for a dataset.  At
#' this stage the summary is prepared for a single specific locus, in contrast
#' to \code{\link{analyze_sample}}.  The Allele1 entries correspond to the
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
#'   \code{\link{analyze_sample}}.
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
  # extract sample data entries that meet all criteria for a potential allele
  locus.name <- sample.attrs[["Locus"]]
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  # Note that counts.locus is more restrictive than the total counts of all
  # entries with MatchingLocus equal to the given locus name, since idx includes
  # the extra checks above.
  count.total <- sum(sample.data$Count)
  count.locus <- sum(chunk$Count)
  # If specified, take the top two matching entries at the given expected
  # lengths.  If there's just one length expected, this should be two entries at
  # that one length.  Otherwise it will be two distinct lengths between the two
  # entries.  We do this in a roundabout way to handle the length-homoplasy
  # case, where only one sequence length is expected but there are two distinct
  # alleles present at that same length.
  expected_lengths <- c(sample.attrs[["ExpectedLength1"]],
                        sample.attrs[["ExpectedLength2"]])
  expected_lengths <- unique(expected_lengths[!is.na(expected_lengths)])
  if (!is.null(expected_lengths) & length(expected_lengths) == 0)
    expected_lengths <- NULL
  if (!is.null(expected_lengths)) {
    idx <- match(expected_lengths, chunk$Length)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 1)
      idx <- c(idx, 1 + match(expected_lengths, chunk$Length[-idx]))
    idx <- idx[!is.na(idx)]
    chunk <- chunk[idx, ]
  }

  # Threshold potential alleles at minimum count.
  # If two expected lengths were given, just take those two entries.
  # If a single expected length was given, ensure we at least get one sequence
  # at that length.
  idx <- chunk$Count >= fraction.min * count.locus
  if (is.null(expected_lengths)) {
    chunk <- chunk[idx, ]
  } else if (length(expected_lengths) == 1) {
  # The only case where we should bother with this, if there are expected
  # lengths given, is if there was only one length.  In that case we should
  # narrow it down to either one or two sequences at that length.
    # (If no rows would be selected, at least take the top entry.)
    if (sum(idx) == 0)
      idx <- 1
    chunk <- chunk[idx, ]
  }

  # Remove ambiguous sequences, if present.
  ambig <- chunk[2, "Ambiguous"]
  chunk <- chunk[! chunk[, "Ambiguous"], ]
  # Remove stutter, if present.
  stutter <- NA
  artifact <- NA
  if (is.null(expected_lengths)) {
    # Remove stutter, if present.
    stutter <- !is.na(chunk[2, "Stutter"])
    chunk <- chunk[is.na(chunk[, "Stutter"]), ]
    # Remove artifact-like sequences, if present.
    artifact <- !is.na(chunk[2, "Artifact"])
    chunk <- chunk[is.na(chunk[, "Artifact"]), ]
  }
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- nrow(chunk)
  # Enforce count limit after all filtering (but before accounting for the
  # stutter removal or fraction thresholding above)
  if (is.null(expected_lengths) && count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[1, attr.names]
  allele2 <- chunk[2, attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  homozygous <- nrow(chunk) == 1
  sample.summary <- c(allele1, allele2,
                      list(Homozygous = homozygous,
                           Ambiguous = ambig,
                           Stutter = stutter,
                           Artifact = artifact,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}


#' Summarize a processed STR sample, Simple Version
#'
#' Summarize as in \code{\link{summarize_sample}}, but skip allele_match and
#' stutter removal.  Entries in the returned list are the same.
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_sample}}.
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
  locus.name <- sample.attrs[["Locus"]]
  chunk <- with(sample.data,
                sample.data[LengthMatch & ! is.na(LengthMatch), ])
  count.total <- sum(sample.data$Count)
  count.locus <- sum(chunk$Count)
  chunk <- chunk[chunk$Count >= fraction.min * count.locus, ]
  prominent.seqs <- nrow(chunk)
  if (count.locus < counts.min) {
    chunk <- chunk[0, ]
  }
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[1, attr.names]
  allele2 <- chunk[2, attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  homozygous <- nrow(chunk) == 1
  sample.summary <- c(allele1, allele2,
                      list(Homozygous = homozygous,
                           Ambiguous = NA,
                           Stutter = NA,
                           Artifact = NA,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

#' Summarize a processed STR sample, Length Version
#'
#' Summarize as in \code{\link{summarize_sample}}, but group entries by sequence
#' length rather than identity.
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_sample}}.
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
summarize_sample_by_length <- function (sample.data, sample.attrs,
                                        fraction.min=0.15,
                                        counts.min=500) {
  locus.name <- sample.attrs[["Locus"]]
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  count.total <- sum(sample.data$Count)
  count.locus <- sum(chunk$Count)
  chunk <- with(chunk, {
         chunk %>%
         dplyr::group_by(Length, MatchingLocus, MotifMatch, LengthMatch) %>%
           dplyr::summarize(Count = sum(Count),
                            Seq = first(Seq),
                            Stutter = first(Stutter))
       })

  chunk <- chunk[order(chunk$Count, decreasing = T), ]
  chunk <- chunk[chunk$Count >= fraction.min * count.locus, ]

  # Remove stutter, if present.
  stutter <- !is.na(chunk[2, "Stutter"])
  chunk <- chunk[is.na(chunk[, "Stutter"]), ]
  # Remove artifact-like sequences, if present.
  artifact <- !is.na(chunk[2, "Artifact"])
  chunk <- chunk[is.na(chunk[, "Artifact"]), ]

  prominent.seqs <- nrow(chunk)
  # enforce count limit after all filtering
  if (count.locus < counts.min)
    chunk <- chunk[0, ]
  # Take top to remaining entries as the two alleles and keep selected
  # attributes.
  attr.names <- c("Seq", "Count", "Length")
  allele1 <- chunk[1, attr.names]
  allele2 <- chunk[2, attr.names]
  colnames(allele1) <- paste0("Allele1", colnames(allele1))
  colnames(allele2) <- paste0("Allele2", colnames(allele2))
  # Combine into summary list with additional attributes.
  homozygous <- nrow(chunk) == 1
  sample.summary <- c(allele1, allele2,
                      list(Homozygous = homozygous,
                           Stutter = stutter,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

#' Check Sample Data for Potential Allele Matches
#'
#' Check the entries in a processed sample data frame for potential matches to a
#' given locus.
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_sample}}.
#' @param locus.name character name of locus to match against.
#'
#' @return logical vector of entries for potential alleles.
allele_match <- function(sample.data, locus.name) {
  with(sample.data,
       as.character(MatchingLocus) == locus.name &
         MotifMatch &
         LengthMatch)
}
