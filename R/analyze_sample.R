#' Analyze sequence table and categorize sequences
#'
#' Converts a full STR sequence data frame into a per-locus version and adds a
#' Category factor column to designate which sequences look like alleles,
#' artifacts, etc.  At this stage the summary is prepared for a single specific
#' locus, in contrast to \code{\link{analyze_seqs}}.
#'
#' @details
#' Factor levels in the added Category column:
#'  * Allele: An identified allele sequence.  There will be between zero and two
#'  of these.
#'  * Prominent: Any additional sequences beyond two called alleles that match
#'  all requirements (sequences that match all locus attributes, do not appear
#'  artifactual, and are above a given fraction of filtered reads).
#'  * Insignificant: Sequences with counts below the \code{fraction.min}
#'  threshold.
#'  * Ambiguous: Sequences passing the \code{fraction.min} threshold but with
#'  non-ACTG characters such as N, as defined by the Ambiguous column.
#'  * Stutter: Sequences passing the \code{fraction.min} threshold but matching
#'  stutter sequence criteria as defined by the Stutter column.
#'  * Artifact: Sequences passing the \code{fraction.min} threshold but matching
#'  non-stutter artifact sequence criteria as defined by the Artifact column.
#' @md
#'
#' @param sample.data data frame of processed data for sample as produced by
#'   \code{\link{analyze_seqs}}.
#' @param sample.attrs list of sample attributes, such as the rows produced by
#'   \code{\link{prepare_dataset}}.  Used to select the locus name to filter on.
#' @param fraction.min numeric threshold for the minimum proportion of counts a
#'   given entry must have, compared to the total matching all criteria for that
#'   locus, to be considered as a potential allele.
#'
#' @return filtered version of \code{sample.data} with added Category column.
analyze_sample <- function(sample.data, sample.attrs, fraction.min) {
  # Extract sample data entries that meet all criteria for a potential allele.
  locus.name <- unlist(sample.attrs["Locus"])
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  within(chunk, {
    Category <- factor(, levels = c("Allele", "Prominent", "Insignificant",
                                    "Ambiguous", "Stutter", "Artifact"))
    # Threshold potential alleles at minimum count
    Category[Count < fraction.min * sum(Count)] <- "Insignificant"
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
  locus.name <- unlist(sample.attrs["Locus"])
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]

  # If specified, take the top two matching entries at the given expected
  # lengths.  If there's just one length expected, this should be two entries at
  # that one length.  Otherwise it will be two distinct lengths between the two
  # entries.  We do this in a roundabout way to handle the length-homoplasy
  # case, where only one sequence length is expected but there are two distinct
  # alleles present at that same length.

  # Tidy up expected lengths.
  expected_lengths <- as.integer(unlist(sample.attrs[c("ExpectedLength1",
                                                       "ExpectedLength2")]))
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
             Category[Count < fraction.min * sum(Count)] <- "Insignificant"
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

analyze_sample_naive <- function(sample.data, sample.attrs, fraction.min) {
  idxl <- with(sample.data, LengthMatch & ! is.na(LengthMatch))
  chunk <- sample.data[idxl, ]
  within(chunk, {
    Category <- factor(, levels = c("Allele", "Prominent", "Insignificant",
                                    "Ambiguous", "Stutter", "Artifact"))
    # Threshold potential alleles at minimum count
    Category[Count < fraction.min * sum(Count)] <- "Insignificant"
    # Top two remaining will be called alleles
    Category[is.na(Category)][0:min(2, sum(is.na(Category)))] <- "Allele"
    # The rest are prominent but not allele-level
    Category[is.na(Category)] <- "Prominent"
  })
}
