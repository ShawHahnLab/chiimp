# Analyze a single set of sequences and return full results.

#' Analyze a set of STR sequences
#'
#' Dereplicates the given sequences and annotates any STR sequences found,
#' returning the processed data as a data frame with one row per unique
#' sequence, sorted by count.  At this stage no information is filtered out, and
#' all loci are treated equally.
#'
#' @details
#' Columns in the returned data frame:
#'  * `Seq`: sequence text for each unique sequence
#'  * `Count`: integer count of occurrences of this exact sequence
#'  * `Length`: integer sequence length
#'  * `MatchingLocus`: factor for the name of the locus matching each
#'    sequence, by checking the primer
#'  * `MotifMatch`: logical: are there are least `nrepeats` perfect
#'    adjacent repeats of the STR motif for the matching locus?
#'  * `LengthMatch`: logical: is the sequence length within the expected
#'    range for the matching locus?
#'  * `Ambiguous`: logical: are there unexpected characters in the sequence
#'    content?
#'  * `Stutter`: integer: for any sequence that looks like potential PCR
#'    stutter, the index of the row that may be the source of the stutter band.
#'  * `Artifact`: integer: for any sequence that looks like potential PCR
#'    artifact (other than stutter), the index of the row that may be the source
#'    of the stutter band.
#'  * `FractionOfTotal`: numeric fraction of the number of sequences
#'    represented by each unique sequence compared to the total.
#'  * `FractionOfLocus`: numeric fraction of the number of sequences
#'    represented by each unique sequence compared to the total for that
#'    particular matching locus.
#'
#' @param seqs character vector containing sequences.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param nrepeats number of repeats of each locus' motif to require for a
#'   match.
#' @param max_stutter_ratio highest ratio of read counts for second most
#'   frequent sequence to the most frequent where the second will be
#'   considered stutter.
#' @param artifact.count.ratio_max as for `max_stutter_ratio` but for
#'   non-stutter artifact sequences.
#' @param ... additional arguments for [make_read_primer_table]
#'
#' @return data frame of dereplicated sequences with added annotations.
#'
#' @export
#'
#' @examples
#' # Starting from non-locus-specific sequences,
#' # a locus attributes table, and requiring
#' # three side-by-side motif repeats to register
#' # as a motif match for a locus,
#' raw_seq_vector <- c(test_data$seqs1$A, test_data$seqs1$B)
#' locus_attrs <- test_data$locus_attrs
#' num_adjacent_repeats <- 3
#' # Convert the character vector of sequences
#' # into a data frame with one row per
#' # unique sequence.
#' seq_data <- analyze_seqs(raw_seq_vector,
#'                          locus_attrs,
#'                          num_adjacent_repeats)
#' @md
analyze_seqs <- function(
  seqs, locus_attrs, nrepeats = cfg("min_motif_repeats"),
  max_stutter_ratio = cfg("max_stutter_ratio"),
  artifact.count.ratio_max = cfg("max_artifact_ratio"),
  ...) {

  primer_info <- make_read_primer_table(seqs, locus_attrs, ...)
  primer_info$MatchingLocus <- factor(
    primer_info$MatchingLocus, levels = unique(locus_attrs$Locus))

  # Dereplicate
  if (nrow(primer_info) > 0) {
    data <- do.call(rbind, lapply(
        split(primer_info, primer_info$Seq), function(chunk) {
      data.frame(
        Seq = chunk$Seq[1],
        Count = nrow(chunk),
        Length = nchar(chunk$Seq[1]),
        MatchingLocus = chunk$MatchingLocus[1],
        stringsAsFactors = FALSE)
    }))
    data <- data[order(data$Count, decreasing = TRUE), ]
  } else {
    data <- data.frame(
      Seq = character(),
      Count = integer(),
      Length = integer(),
      MatchingLocus = factor(c(), levels = unique(locus_attrs$Locus)),
      stringsAsFactors = FALSE)
  }
  rownames(data) <- NULL

  # Label rows where the sequence content matches the matching locus" repeat
  # motif.
  data$MotifMatch <- check_motif(data, locus_attrs, nrepeats)
  # Label rows where the sequence length matches the matching locus" length
  # range.
  data$LengthMatch <- check_length(data, locus_attrs)
  # Label rows that contain unexpected characters in the sequence content.
  data$Ambiguous <- ! grepl("^[ACGT]*$", data$Seq, ignore.case = TRUE)
  # Label rows that look like PCR stutter or other artifacts of other rows.
  data$Stutter <- find_stutter(data, locus_attrs, max_stutter_ratio)
  data$Artifact <- find_artifact(data, locus_attrs, artifact.count.ratio_max)
  # Add columns for the proportion of counts out of the total and out of those
  # for the matching locus.  This way this information is preserved even in
  # subsets of the original sample data.
  data$FractionOfTotal <- data$Count / sum(data$Count)
  data$FractionOfLocus <- with(data, {
    total_per_locus <- unlist(lapply(levels(MatchingLocus), function(loc) {
      sum(data[MatchingLocus == loc, "Count"], na.rm = TRUE)
      }))
    names(total_per_locus) <- levels(MatchingLocus)
    Count / total_per_locus[MatchingLocus]
  })
  return(data)
}

#' Check sequences for STR repeats
#'
#' Return a logical vector specifying, for each entry, if the sequence contains
#' at least `nrepeats` perfect repeats of the matching locus' motif.
#'
#' @param sample.data data frame of processed sequence data.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param nrepeats number of repeats of each locus" motif to require for a
#'   match.
#'
#' @return logical vector indicating where repeats were observed.
#' @md
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
#' @param sample.data data frame of processed sequence data.
#' @param locus_attrs data frame of attributes for loci to look for.
#'
#' @return logical vector specifying, for each entry, if the sequence is within
#'   the matching locus' expected length range.
#' @md
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

#' Check for PCR stutter between sequences
#'
#' Searches a processed STR sample for entries that may be artifacts of PCR
#' stutter from another entry in the sample.  This only considers STR-labeled
#' rows and requires a given entry to have counts at most `count.ratio_max`
#' compared to the candidate "source" entry to be considered stutter.  Sequence
#' content is not currently considered, just relative sequence lengths and
#' counts.
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param count.ratio_max comparing the currently-checked entry to another
#'   entry, this is the highest ratio of counts where an entry will still be
#'   considered stutter.
#'
#' @return integer vector specifying, for each entry, the row index for another
#'   entry that may have produced each entry as a stutter band.
#' @md
find_stutter <- function(
  sample.data, locus_attrs, count.ratio_max = cfg("max_stutter_ratio")) {

  stutter <- integer(nrow(sample.data)) * NA

  # across rows for this locus,
  for (locus_name in rownames(locus_attrs)) {
    idxl_main <- sample.data$MatchingLocus == locus_name
    idxl_main[is.na(idxl_main)] <- FALSE
    d <- sample.data[idxl_main, ]
    motif.len <- nchar(as.character(locus_attrs[locus_name, "Motif"]))
    # Any given index in d could be stutter.
    # these are lengths to look for as potential stutter sources of each row
    lens_up <- d$Length + motif.len
    # These will be row indexes for any possible stutter source, going by
    # length.  Implicitly each index will match the highest count as well.
    idx_length <- match(lens_up, d$Length)
    # TRUE for the idx_length rows with counts below the cutoff
    idxl_stutter <- d$Count < d$Count[idx_length] * count.ratio_max
    # Now, these are the index values in d that look like stutter of any given
    # row (non-matching entries are left in place but set to NA)
    idx_length[! idxl_stutter] <- NA
    # Next, insert these into the stutter vector, at the right position, with
    # the right index values.
    stutter[which(idxl_main)] <- idx_length + cumsum(!idxl_main)[idx_length]
  }

  return(stutter)
}

#' Check for non-stutter PCR artifacts between sequences
#'
#' Searches a processed STR sample for entries that may be PCR artifacts, other
#' than stutter, from another entry in the sample.  Potential artifacts are
#' sequences with counts lower than another sequence by a given ratio and
#' sequence length within 1 nucleotide of the other sequence.  This only
#' considers STR-labeled rows and requires a given entry to have counts at most
#' `count.ratio_max` compared to the candidate "source" entry to be considered
#' an artifact.  Sequence content is not currently considered, just relative
#' sequence lengths and counts.
#'
#' @param sample.data data frame of processed sample data.
#' @param locus_attrs data frame of attributes for loci to look for.
#' @param count.ratio_max comparing the currently-checked entry to another
#'   entry, this is the highest ratio of counts where an entry will still be
#'   considered artifactual
#'
#' @return integer vector specifying, for each entry, the row index for another
#'   entry that may have produced each entry as an artifactual sequence.
#' @md
find_artifact <- function(
    sample.data, locus_attrs, count.ratio_max = cfg("max_artifact_ratio")) {

  artifact <- integer(nrow(sample.data)) * NA

  # across rows for this locus,
  for (locus_name in rownames(locus_attrs)) {
    idxl_main <- sample.data$MatchingLocus == locus_name
    idxl_main[is.na(idxl_main)] <- FALSE
    d <- sample.data[idxl_main, ]
    # for each length it's the index of the first (highest-count) entry with
    # that length.
    matches <- matrix(c(match(d$Length, d$Length),
                        match(d$Length - 1, d$Length),
                        match(d$Length + 1, d$Length)), ncol = 3)
    # Now take the first (highest-count) of the set, ignoring NA's, and giving
    # NA if all is NA for a given entry.
    idx_length <- apply(matches, 1, function(r) {
      if (all(is.na(r))) {
        return(NA)
      }
      min(r, na.rm = TRUE)
    })
    # TRUE for the idx_length rows with counts below the cutoff
    idxl_artifact <- d$Count < d$Count[idx_length] * count.ratio_max
    # Now, these are the index values in d that look like stutter of any given
    # row (non-matching entries are left in place but set to NA)
    idx_length[! idxl_artifact] <- NA
    # Next, insert these into the stutter vector, at the right position, with
    # the right index values.
    artifact[which(idxl_main)] <- idx_length + cumsum(!idxl_main)[idx_length]
  }

  return(artifact)
}
