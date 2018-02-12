# Analyze a single sample and return full results.

#' Analyze a set of STR sequences
#'
#' Dereplicates the given sequences and annotates any STR sequences found,
#' returning the processed sample as a data frame with one row per unique
#' sequence, sorted by count.  At this stage no information is filtered out, and
#' all loci are treated equally.
#'
#' @details
#' Columns in the returned data frame:
#'  * Seq: sequence text for each unique sequence
#'  * Count: integer count of occurrences of this exact sequence
#'  * Length: integer sequence length
#'  * MatchingLocus: factor for the name of the locus matching each sequence,
#'    by checking the primer
#'  * MotifMatch: logical -- are there are least \code{nrepeats} perfect
#'    adjacent repeats of the STR motif for the matching locus?
#'  * LengthMatch: logical -- is the sequence length within the expected range
#'    for the matching locus?
#'  * Stutter: integer: for any sequence that looks like potential PCR stutter,
#'    the index of the row that may be the source of the stutter band.
#'  * FractionOfTotal: numeric fraction of the number of sequences
#'    represented by each unique sequence compared to the total.
#'  * FractionOfLocus: numeric fraction of the number of sequences represented
#'    by each unique sequence compared to the total for that particular
#'    matching locus.
#' @md
#'
#' @param seqs character vector containing sequences.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match.
#'
#' @return data frame of dereplicated sequences with added annotations.
#'
#' @export
analyze_sample <- function(seqs, locus_attrs, nrepeats) {
  # Dereplicate sequences
  tbl <- table(seqs)
  count <- as.integer(tbl)
  seqs <- as.character(names(tbl))
  # Combine into initial data frame.  Leave character vectors as-is here since
  # they"re unique anyway.  Sort by number of occurrences, decreasing. Renumber
  # rows.
  data <- data.frame(Seq = seqs,
                     Count = count,
                     Length = nchar(seqs),
                     stringsAsFactors = FALSE)
  data <- data[order(data$Count, decreasing = T), ]
  rownames(data) <- NULL
  # Label rows with the apparent locus by checking primer sequences.  Note that
  # this uses the first matching locus for each row.
  data$MatchingLocus <- find_matching_primer(data, locus_attrs)
  # Label rows where the sequence content matches the matching locus" repeat
  # motif.
  data$MotifMatch <- check_motif(data, locus_attrs, nrepeats)
  # Label rows where the sequence length matches the matching locus" length
  # range.
  data$LengthMatch <- check_length(data, locus_attrs)
  # Label rows that look like PCR stutter of other rows.
  data$Stutter <- find_stutter(data, locus_attrs)
  # Add columns for the proportion of counts out of the total and out of those
  # for the matching locus.  This way this information is preserved even in
  # subsets of the original sample data.
  data$FractionOfTotal <- data$Count / sum(data$Count)
  data$FractionOfLocus <- with(data, {
    total_per_locus <- unlist(lapply(levels(MatchingLocus), function(loc)
      sum(data[MatchingLocus == loc, "Count"], na.rm = T)))
    names(total_per_locus) <- levels(MatchingLocus)
    Count / total_per_locus[MatchingLocus]
  })
  return(data)
}

#' Label STR sample rows by locus
#'
#' Return a factor giving the locus name matching each sample data entry's
#' matching primer (first found).
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#'
#' @return factor of locus names corresponding the matched primer sequences.
find_matching_primer <- function(sample.data, locus_attrs) {
  # Separately check each primer.  Is there a slicker way to do this all at
  # once?
  matches <- do.call(cbind, lapply(rownames(locus_attrs), function(locus_name) {
    primer <- as.character(locus_attrs[locus_name, "Primer"])
    result <- grepl(primer, substr(sample.data$Seq, 1, nchar(primer) + 10))
    c(locus_name)[as.numeric((!result) + 1)]
  }))
  # Collapse that set down to just the first match for each entry.
  first.matches <- apply(matches, 1, function(m) m[match(TRUE, !is.na(m))])
  # Store that as a factor and use the standard locus name levels.
  factor(first.matches, levels = rownames(locus_attrs))
}

#' Check sequences for STR repeats
#'
#' Return a logical vector specifying, for each entry, if the sequence contains
#' at least \code{nrepeats} perfect repeats of the matching locus' motif.
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param nrepeats number of repeats of each locus" motif to require for a
#'   match.
#'
#' @return logical vector indicating where repeats were observed.
check_motif <- function(sample.data, locus_attrs, nrepeats) {
  with(sample.data, {
    motif <- locus_attrs[MatchingLocus, "Motif"]
    # this is kind of awful.  clean this up somehow.
    locus_repeat <- unlist(lapply(data.frame(matrix(
      rep(motif, each = nrepeats), nrow = nrepeats)), paste, collapse = ""))
    pairs <- cbind(locus_repeat, Seq)
    apply(pairs, 1, function(p) grepl(p[1], p[2]))
  })
}

#' Check sequences for expected length
#'
#' Return a logical vector specifying, for each entry, if the sequence is within
#' the expected length range for the matched locus.
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#'
#' @return logical vector specifying, for each entry, if the sequence is within
#'   the matching locus' expected length range.
check_length <- function(sample.data, locus_attrs) {
  with(sample.data, {
    Lmin <- locus_attrs[MatchingLocus, "LengthMin"]
    Lmax <- locus_attrs[MatchingLocus, "LengthMax"]
    Lbuff <- locus_attrs[MatchingLocus, "LengthBuffer"]
    lmin <- Lmin - Lbuff
    lmax <- Lmax + Lbuff
    Length >= lmin & Length <= lmax
  })
}

#' Check for PCR stutter between sample entries
#'
#' Searches a processed STR sample for entries that may be artifacts of PCR
#' stutter from another entry in the sample.  This only considers STR-labeled
#' rows and those with a count above \code{count.min}, and requires a given
#' entry to have counts at most \code{count.ratio_max} compared to the candidate
#' "source" entry to be considered stutter.  Sequence content is not currently
#' considered, just relative sequence lengths and counts.
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param count.ratio_max comparing the currently-checked entry to another
#'   entry, this is the highest ratio of counts where an entry will still be
#'   considered stutter.
#' @param count.min lowest count for a given entry to still attempt a check for
#'   stutter.
#'
#' @return integer vector specifying, for each entry, the row index for another
#'   entry that may have produced each entry as a stutter band.
find_stutter <- function(sample.data, locus_attrs,
                         count.ratio_max = 1 / 3,
                         count.min = 10) {

  stutter <- as.integer(rep(NA, nrow(sample.data)))

  # across rows for this locus,
  for (locus_name in rownames(locus_attrs)) {
    locus.match <- sample.data$MatchingLocus == locus_name
    motif.len <- nchar(as.character(locus_attrs[locus_name, "Motif"]))
    for (idx in which(locus.match)) {
      # look at each row, and see if any others may be
      # stutter using each row here as a source.
      L <- sample.data[idx, "Length"] - motif.len
      C <- sample.data[idx, "Count"] * count.ratio_max
      #  mark any of the other rows that look like stutter of this row.
      others.match <- locus.match &
        sample.data$Length == L &
        sample.data$Count <= C &
        sample.data$Count > count.min
      stutter[others.match] <- idx
    }
  }

  return(stutter)
}
