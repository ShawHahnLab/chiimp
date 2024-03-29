#' Versions of analyze_sample
#'
#' This list shows the different versions of [analyze_sample]
#' recognized for use via configuration in full_analysis.  See also
#' [sample_summary_funcs] and [CFG_DEFAULTS].
#'
#' @export
#' @md
sample_analysis_funcs <- c("analyze_sample",
                           "analyze_sample_guided",
                           "analyze_sample_naive")

#' Analyze sequence table and categorize sequences
#'
#' Converts a full STR sequence data frame into a per-locus version and adds a
#' Category factor column to designate which sequences look like alleles,
#' artifacts, etc.  At this stage the summary is prepared for a single specific
#' locus, in contrast to [analyze_seqs].  See the Details section
#' below for a description of the factor levels in the new Category column, and
#' see the Functions section below for how specific variants of this function
#' behave.
#'
#' @details
#' Factor levels in the added Category column, in order:
#'  * Allele: An identified allele sequence.  There will be between zero and two
#'  of these.
#'  * Prominent: Any additional sequences beyond two called alleles that match
#'  all requirements (sequences that match all locus attributes, do not appear
#'  artifactual, and are above a given fraction of filtered reads).
#'  * Insignificant: Sequences with counts below the `min_allele_abundance`
#'  threshold.
#'  * Ambiguous: Sequences passing the `min_allele_abundance` threshold but with
#'  non-ACTG characters such as N, as defined by the Ambiguous column of
#'  `seq_data`.
#'  * Stutter: Sequences passing the `min_allele_abundance` threshold but
#'  matching stutter sequence criteria as defined by the Stutter column of
#'  `seq_data`.
#'  * Artifact: Sequences passing the `min_allele_abundance` threshold but
#'  matching non-stutter artifact sequence criteria as defined by the Artifact
#'  column of `seq_data`.
#' @md
#'
#' @param seq_data data frame of processed data for sample as produced by
#'   [analyze_seqs].
#' @param sample_attrs list of sample attributes, such as the rows produced by
#'   [prepare_dataset].  Used to select the locus name to filter on.
#' @param min_allele_abundance numeric threshold for the minimum proportion of
#'   counts a given entry must have, compared to the total matching all criteria
#'   for that locus, to be considered as a potential allele.
#'
#' @return filtered version of `seq_data` with added Category column.
#'
#' @describeIn analyze_sample default version of sample analysis.  From here use
#'   [summarize_sample].
#'
#' @export
#' @md
analyze_sample <- function(
  seq_data, sample_attrs, min_allele_abundance = cfg("min_allele_abundance")) {
  # Extract sample data entries that meet all criteria for a potential allele.
  locus_name <- unlist(sample_attrs["Locus"])
  idx <- which(allele_match(seq_data, locus_name))
  chunk <- seq_data[idx, ]
  attr(chunk, "min_allele_abundance") <- min_allele_abundance
  within(chunk, {
    Category <- factor(, levels = c("Allele", "Prominent", "Insignificant",
                                    "Ambiguous", "Stutter", "Artifact"))
    # Threshold potential alleles at minimum count
    Category[Count < min_allele_abundance * sum(Count)] <- "Insignificant"
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

#' @describeIn analyze_sample version of sample analysis guided by expected
#'   sequence length values.  Additional items `ExpectedLength1` and
#'   optionally `ExpectedLength2` can be supplied in the
#'   `sample_attrs` list.  If NA or missing the behavior will match
#'   [analyze_sample].  If two expected lengths are given, the
#'   min_allele_abundance argument is ignored.  If at least one expected length
#'   is given, the stutter/artifact filtering is disabled.  From here use
#'   [summarize_sample_guided].
#'
#' @export
#' @md
analyze_sample_guided <- function(
  seq_data, sample_attrs, min_allele_abundance = cfg("min_allele_abundance")) {
  # Extract sample data entries that meet all criteria for a potential allele.
  locus_name <- unlist(sample_attrs["Locus"])
  idx <- which(allele_match(seq_data, locus_name))
  chunk <- seq_data[idx, ]
  attr(chunk, "min_allele_abundance") <- min_allele_abundance

  # If specified, take the top two matching entries at the given expected
  # lengths.  If there's just one length expected, this should be two entries at
  # that one length.  Otherwise it will be two distinct lengths between the two
  # entries.  We do this in a roundabout way to handle the length-homoplasy
  # case, where only one sequence length is expected but there are two distinct
  # alleles present at that same length.

  # Tidy up expected lengths.
  expected_lengths <- as.integer(unlist(sample_attrs[c("ExpectedLength1",
                                                       "ExpectedLength2")]))
  expected_lengths <- unique(expected_lengths[! is.na(expected_lengths)])

  categories <- c("Allele", "Prominent", "Insignificant",
                  "Ambiguous", "Stutter", "Artifact")

  switch(length(expected_lengths) + 1,
         # Zero expected lengths: analyze as usual
         analyze_sample(seq_data, sample_attrs, min_allele_abundance), {
         # One expected length: may be homozygous or heterozygous.
           # Find rows of interest, matching expected length.
           idxl <- chunk$Length == expected_lengths
           within(chunk, {
             Category <- factor(, levels = categories)
             # Exclude ambiguous sequences first.
             Category[Ambiguous] <- "Ambiguous"
             # Assign top seq at expected length as first allele.
             Category[is.na(Category) & idxl][1] <- "Allele"
             # Threshold potential alleles at minimum count.
             Category[
               is.na(Category) & Count < min_allele_abundance * sum(Count)] <-
               "Insignificant"
             # Assign second seq at expected length as second allele, if there
             # is one.
             Category[is.na(Category) & idxl][1] <- "Allele"
             # And that's it.  We make no comment on the remaining entries and
             # leave them as NA.
           })
         }, {
         # Two expected lengths: definitely heterozygous.  No need to consider
         # fractions here.
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

#' @describeIn analyze_sample version of sample analysis without
#'   stutter/artifact filtering.  From here use [summarize_sample] as for
#'   `analyze_sample`.
#'
#' @export
#' @md
analyze_sample_naive <- function(
  seq_data, sample_attrs, min_allele_abundance = cfg("min_allele_abundance")) {
  idxl <- with(seq_data, LengthMatch & ! is.na(LengthMatch))
  chunk <- seq_data[idxl, ]
  attr(chunk, "min_allele_abundance") <- min_allele_abundance
  within(chunk, {
    Category <- factor(, levels = c("Allele", "Prominent", "Insignificant",
                                    "Ambiguous", "Stutter", "Artifact"))
    # Threshold potential alleles at minimum count
    Category[Count < min_allele_abundance * sum(Count)] <- "Insignificant"
    # Top two remaining will be called alleles
    Category[is.na(Category)][0:min(2, sum(is.na(Category)))] <- "Allele"
    # The rest are prominent but not allele-level
    Category[is.na(Category)] <- "Prominent"
  })
}

#' Check Sample Data for Potential Allele Matches
#'
#' Check the entries in a processed sample data frame for potential matches to a
#' given locus.
#'
#' @param sample_data data frame of processed data for sample as produced by
#'   [analyze_seqs].
#' @param locus_name character name of locus to match against.
#'
#' @return logical vector of entries for potential alleles.
#' @md
allele_match <- function(sample_data, locus_name) {
  with(sample_data,
       as.character(MatchingLocus) == locus_name &
         MotifMatch &
         LengthMatch)
}
