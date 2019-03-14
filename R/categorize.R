# Interpret genotyping results for samples with known identity.

#' Associate known genotypes with samples
#'
#' Using the Name column of the given results summary data frame, pair each
#' called genotype with the known alleles.  A data frame with two columns,
#' \code{CorrectAllele1Seq} and \code{CorrectAllele2Seq}, is returned. If
#' matching entries are found in \code{Allele1Seq} and/or \code{Allele2Seq} the
#' order will be preserved, and at this point the two allele entries should
#' match up directly for genotypes that were called correctly.
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param genotypes.known data frame of known genotypes that should be compared
#'   to the observed genotypes in the results, as loaded by
#'   \code{\link{load_genotypes}}.
#'
#' @return data frame with two columns for the two correct alleles, and rows
#'   matching the input summary table.
#'
#' @export
match_known_genotypes <- function(results_summary, genotypes.known) {
  # match name/locus combos with genotypes
  id_tbl <- paste(results_summary$Name, results_summary$Locus)
  id_kg <- paste(genotypes.known$Name, genotypes.known$Locus)
  idx <- match(id_tbl, id_kg)
  # Build data frame of correct allele sequences
  result <- data.frame(CorrectAllele1Seq = genotypes.known[idx, "Allele1Seq"],
                       CorrectAllele2Seq = genotypes.known[idx, "Allele2Seq"],
                       stringsAsFactors = FALSE)
  # Ensure ordering within pairs matches samples, if possible.
  for (i in 1:nrow(result)) {
    a <- results_summary[i, c("Allele1Seq", "Allele2Seq")]
    kg <- result[i, ]
    idx <- match(a, kg)
    if (idx[1] %in% 2 || idx[2] %in% 1)
      result[i, ] <- rev(kg)
  }
  result
}

#' Categorize genotyping results
#'
#' For a given results summary data frame that has \code{CorrectAllele1Seq} and
#' \code{CorrectAllele2Seq} columns (such as produced by
#' \code{\link{match_known_genotypes}}) added, create a factor labeling every
#' row of the input data frame by its genotyping outcome.
#'
#' @details
#' Levels in the returned factor, in order:
#'
#' * Correct: one/two alleles match.
#' * Incorrect at least one allele does not match.
#' * Blank: No alleles were called in the analysis even though known genotypes
#'    were supplied.
#' * Dropped Allele: One called allele is correct for a heterozygous individual,
#'   but no second allele was called.
#'
#' Cases that should not occur, such as \code{CorrectAllele1Seq} and
#' \code{CorrectAllele2Seq} both set to NA, map to NA in the returned factor.
#' @md
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}} with extra columns as produced by
#'   \code{\link{match_known_genotypes}}.
#'
#' @return factor defining genotyping result category for every row of the input
#'   data frame.
#'
#' @export
categorize_genotype_results <- function(results_summary) {
  # Five possibilities for either NA/not NA plus outcome of non-NA pair
  # All five possibilities for a single allele check:
  #   0: Both non-NA, simple mismatch
  #   1: A not NA, C NA (no correct allele matched this one)
  #   2: A NA, C not NA (we missed a correct allele and left this blank)
  #   3: A NA, C NA (correctly did not report an allele)
  #   4: Both non-NA, match
  check_allele <- function(allele, ref) {
    a <- is.na(allele) * 2 + is.na(ref) # NA: 1, not NA: 0
    a[a == 0 & allele == ref] <- 4 # special distinction for one case
    a
  }

  # Now, combine for both alleles to have all possible outcomes, and offset by
  # one to account for R's indexing.
  a1 <- check_allele(results_summary$Allele1Seq,
                     results_summary$CorrectAllele1Seq)
  a2 <- check_allele(results_summary$Allele2Seq,
                     results_summary$CorrectAllele2Seq)
  a <- a1 * 5 + a2 + 1

  # Here's all the possible outcomes, categorized.  Cases that should never come
  # up for correctly-labeled genotypes will evaluate to NA.
  lvls <- c(
    # A1 0: first allele simple mismatch.  Whatever A2 is, this is Incorrect.
    "Incorrect", # both mismatch
    "Incorrect", # extra allele, mismatch
    "Incorrect", # drop
    "Incorrect", # correctly missing
    "Incorrect", # second correct
    # A1 1: first allele called, but no correct allele listed.  Still Incorrect.
    "Incorrect", # simple mismatch
    NA,          # second allele also not present?? weird case
    "Incorrect", # both mismatch
    NA,          # no correct allele listed for second either?? weird case
    "Incorrect", # second is correct but first was wrong
    # A1 2: first allele incorrectly blank.
    "Incorrect", # simple mismatch
    "Incorrect", # wrong
    "Blank",     # second allele also incorrectly blank
    "Incorrect", # though this *was* homozygous; we at least got that right.
    "Dropped Allele", # Got one right, but missed A1.
    # A1 3: first allele correctly blank (expecting true homozygote).
    "Incorrect", # simple mismatch
    NA,          # but C2 also NA? weird case
    "Blank",     # A2 also blank
    NA,          # A2 NA but C2 also NA? weird case
    "Correct",   # correct homozygote
    # A1 4: first allele correct.
    "Incorrect", # but second wrong.
    "Incorrect", # second wrongly given when should be blank.
    "Dropped Allele", # Got one right, but missed A2.
    "Correct",   # correctly did not report a second allele (homozygote)
    "Correct"    # correctly did report a second allele (heterozygote)
  )

  # Map the integers for each case to text categories and create factor.
  factor(lvls[a], levels = c("Correct", "Dropped Allele", "Blank", "Incorrect"))
}
