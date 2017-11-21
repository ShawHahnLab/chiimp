# Create a summary list for a single sample, targeting a single locus.  See
# summarize_sample below for the entries in the returned list.

# Different versions summarize_sample to allow via configuration in
# full_analysis.
sample_summary_funcs <- c("summarize_sample",
                          "summarize_sample_naive",
                          "summarize_sample_by_length")

#' Summarize a processed STR sample
#'
#' Converts a full STR sample data frame into a concise list of consistent
#' attributes, suitable for binding together across samples for a dataset.  At
#' this stage the summary is prepared for a single specific locus, in contrast
#' to \code{analyze_sample}.  The Allele1 entries correspond to the sequence
#' with the highest count, Allele2 the second highest.
#'
#' @details
#' Entries in the returned list:
#'  * For Allele1 and Allele2:
#'    * Seq: sequence text for each allele.
#'    * Count: integer count of occrrences of this exact sequence.
#'    * Length: integer sequence length.
#'  * Homozygous: If the sample appears homozygous (if so, the Allele2 entries
#'  will be NA).
#'  * Stutter: If a potential allele was ignored due to apparent PCR stutter.
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
#'   \code{analyze_sample}
#' @param locus.name character name of locus to summarize with
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
summarize_sample <- function(sample.data, locus.name, fraction.min,
                             counts.min) {
  # extract sample data entries that meet all criteria for a potential allele
  idx <- which(allele_match(sample.data, locus.name))
  chunk <- sample.data[idx, ]
  # Note that counts.locus is more restrictive than the total counts of all
  # entries with MatchingLocus equal to the given locus name, since idx includes
  # the extra checks above.
  count.total <- sum(sample.data$Count)
  count.locus <- sum(chunk$Count)
  # Threshold potential alleles at minimum count
  chunk <- chunk[chunk$Count >= fraction.min * count.locus, ]
  # Remove stutter, if present.
  stutter <- FALSE
  if (!is.na(chunk[2, "Stutter"])) {
    chunk <- chunk[-2, ]
    stutter <- TRUE
  }
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
                           Stutter = stutter,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

# Like summarize_sample above, but skips allele_match and stutter removal.
summarize_sample_naive <- function(sample.data, locus.name, fraction.min,
                                   counts.min) {
  chunk <- with(sample.data,
                sample.data[!is.na(MatchingLocus) &
                              as.character(MatchingLocus) == locus.name, ])
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
                           Stutter = NA,
                           CountTotal = count.total,
                           CountLocus = count.locus,
                           ProminentSeqs = prominent.seqs))
  return(sample.summary)
}

# old version of summary algorithm
summarize_sample_by_length <- function (sample.data, locus.name,
                                        fraction.min=0.15,
                                        counts.min=500) {
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
  stutter <- FALSE
  if (!is.na(chunk[2, "Stutter"])) {
    chunk <- chunk[-2, ]
    stutter <- TRUE
  }
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

allele_match <- function(sample.data, locus.name) {
  with(sample.data,
       as.character(MatchingLocus) == locus.name &
         MotifMatch &
         LengthMatch)
}
