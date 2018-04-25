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
  sample.summary <- c(allele1, allele2,
                      list(Homozygous    = sum(chunk$Category == "Allele") == 1,
                           Ambiguous     = "Ambiguous" %in% chunk$Category[2],
                           Stutter       = "Stutter"   %in% chunk$Category[2],
                           Artifact      = "Artifact"  %in% chunk$Category[2],
                           CountTotal    = sum(sample.data$Count),
                           CountLocus    = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

# TODO finish new version here
# Take a data frame produced by analyze_seqs, restrict to a given locus, and add
# a Category column to classify individual sequences.
analyze_sample <- function(sample.data, sample.attrs, fraction.min) {
  # Extract sample data entries that meet all criteria for a potential allele.
  locus.name <- sample.attrs[["Locus"]]
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  # Note that counts.locus is more restrictive than the total counts of all
  # entries with MatchingLocus equal to the given locus name, since idx includes
  # the extra checks above.
  count.locus <- sum(chunk$Count)
  within(chunk, {
    Category <- factor(, levels = c("Allele", "Prominent", "Insignificant",
                                    "Ambiguous", "Stutter", "Artifact"))
    # Threshold potential alleles at minimum count
    Category[Count < fraction.min * count.locus] <- "Insignificant"
    # Label ambiguous, stutter, artifact sequences, if present.
    Category[is.na(Category) & Ambiguous] <- "Ambiguous"
    Category[is.na(Category) & ! is.na(Stutter)] <- "Stutter"
    Category[is.na(Category) & ! is.na(Artifact)] <- "Artifact"
    # Top two remaining will be called alleles
    Category[is.na(Category)][0:min(2, sum(is.na(Category)))] <- "Allele"
    # The rest are prominent but not allele-level
    Category[is.na(Category)] <- "Prominent"
  })
}

analyze_sample_guided <- function(sample.data, sample.attrs, fraction.min) {
  # Extract sample data entries that meet all criteria for a potential allele.
  locus.name <- sample.attrs[["Locus"]]
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]

  # If specified, take the top two matching entries at the given expected
  # lengths.  If there's just one length expected, this should be two entries at
  # that one length.  Otherwise it will be two distinct lengths between the two
  # entries.  We do this in a roundabout way to handle the length-homoplasy
  # case, where only one sequence length is expected but there are two distinct
  # alleles present at that same length.

  # Tidy up expected lengths.
  expected_lengths <- c(sample.attrs[["ExpectedLength1"]],
                        sample.attrs[["ExpectedLength2"]])
  expected_lengths <- unique(expected_lengths[! is.na(expected_lengths)])

  categories <- c("Allele", "Prominent", "Insignificant",
                  "Ambiguous", "Stutter", "Artifact")

  switch(length(expected_lengths) + 1,
         # Zero expected lengths: analyze as usual
         analyze_sample(sample.data, sample.attrs, fraction.min),
         # One expected length: may be homozygous or heterozygous.
         {
           # Find rows of interest, matching expected length.
           idxl <- chunk$Length == expected_lengths
           within(chunk, {
             Category <- factor(, levels = categories)
             # Exclude ambiguous sequences first.
             Category[Ambiguous] <- "Ambiguous"
             # Assign top seq at expected length as first allele.
             Category[idxl & is.na(Category)][1] <- "Allele"
             # Threshold potential alleles at minimum count.
             Category[Count < fraction.min * count.locus] <- "Insignificant"
             # Assign second seq at expected length as second allele, if there
             # is one.
             Category[idxl & is.na(Category)][1] <- "Allele"
             # And that's it.  We make no comment on the remaining entries and
             # leave them as NA.
           })
         },
         # Two expected lengths: definitely heterozygous.  No need to consider
         # fractions here.
         {
           within(chunk, {
             Category <- factor(, levels = categories)
             # Exclude ambiguous sequences first.
             Category[Ambiguous] <- "Ambiguous"
             # The highest-count non-ambiguous sequences at the expected lengths
             # will be marked as the two alleles.
             Category[! Ambiguous &
                        Length == expected_lengths[1]
                      ][1] <- "Allele"
             Category[! Ambiguous &
                        Length == expected_lengths[2]
                      ][1] <- "Allele"
             # And that's it.  We make no comment on the remaining entries and
             # leave them as NA.
           })
         }
         )
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
  expected_lengths <- c(sample.attrs[["ExpectedLength1"]],
                        sample.attrs[["ExpectedLength2"]])
  expected_lengths <- unique(expected_lengths[! is.na(expected_lengths)])
  
  chunk <- analyze_sample_guided(sample.data, sample.attrs, fraction.min)
  
  # How many entries ended up above the threshold, after all filtering?  Ideally
  # this will be either one or two.
  prominent.seqs <- sum(chunk$Category %in% c("Allele", "Prominent"))
  # Enforce count limit after all filtering
  count.locus <- sum(chunk$Count)
  if (is.null(expected_lengths) && count.locus < counts.min) {
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
  sample.summary <- c(allele1, allele2,
                      list(Homozygous    = sum(chunk$Category == "Allele") == 1,
                           Ambiguous     = "Ambiguous" %in% chunk$Category[2],
                           Stutter       = "Stutter"   %in% chunk$Category[2],
                           Artifact      = "Artifact"  %in% chunk$Category[2],
                           CountTotal    = sum(sample.data$Count),
                           CountLocus    = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}


#' Summarize a processed STR sample, Simple Version
#'
#' Summarize as in \code{\link{summarize_sample}}, but skip allele_match and
#' stutter removal.  Entries in the returned list are the same.
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
