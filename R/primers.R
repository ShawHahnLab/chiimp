# Match Primers -----------------------------------------------------------


#' Make data frame of primer info for reads
#' 
#' `make_read_primer_table` take read sequences and per-locus primer
#' information and produces a data frame of matching information for each read.
#' Mismatches are allowed, but not indels.  The rows in the output data frame
#' will correspond exactly to the read vector given.
#' 
#' Output Columns:
#'
#'  * SeqOrig: Original read sequence
#'  * FwdStart: start position of forward primer match
#'  * Fwdstop: end position of forward primer match
#'  * FwdMismatches: number of mismatches for the forward primer
#'  * FwdLocus: name of the Locus whose forward primer matched
#'  * RevStart: start position of reverse primer match
#'  * RevStop: end position of reverse primer match
#'  * RevMismatches: number of mismatches in the reverse primer
#'  * RevLocus: name of the Locus whose reverse primer matched
#'  * MatchingLocus: name of Locus matching in forward or both (if
#'    `use_reverse_primers`) directions
#'  * Seq: modified sequence based on `primer_action` argument(s), if
#'    applicable
#' 
#' Optionally reads can be modified based on the matched primer sequences in one
#' or both directions; see the `primer_action` arguments for those options.  All
#' positions in the returned data frame are indexed from 1 and are oriented in
#' the same direction as the unmodified reads.
#' 
#' The comparison of reads to primers is not exhaustive.  This method gathers
#' exact matches first and searches the remaining reads for inexact matches in
#' order of decreasing abundance of exact matches, taking the first
#' sufficiently-similar match for each read.  For `max_mismatches` values well
#' below the similarity between any two primer sequences (as will typically be
#' the case) this should result in the best-matching primer being reported for
#' each read, but is much faster than doing an exhaustive all-to-all comparison.
#' 
#' @param seqs character vector of read sequences
#' @param locus_attrs data frame of locus attributes
#' @param max_mismatches integer number of mismatches allowed when checking
#'   primers against reads
#' @param primer_action how should reads be modified based on matched primers?
#'   Can be "none", "keep" (trim sequences to the primers but keep the primers
#'   as-is), "replace" (trim sequences *and* replace the matched primer region
#'   with the sequence from `locus_attrs`), or "remove" (trim sequences to
#'   exclude the primer region).
#' @param primer_action_fwd a `primer_action` value for the forward primer
#' @param primer_action_rev a `primer_action` value for the reverse primer
#' @param max_mismatches_fwd a `max_mismatches` value for the forward primer
#' @param max_mismatches_rev a `max_mismatches` value for the reverse primer
#' @param use_reverse_primers use reverse primers, or just forward?
#' @param reverse_primer_r1 is each reverse primer in `locus_attrs` given in its
#'   orientation on the forward read?  This is used to determine how the primers
#'   and reads should be reverse complemented before comparing.
#' @returns data frame of matched primer information, one row for each input
#'   sequence
#' @export
#' @md
make_read_primer_table <- function(
  seqs, locus_attrs,
  max_mismatches=config.defaults$seq_analysis$max_mismatches,
  primer_action=config.defaults$seq_analysis$primer_action,
  primer_action_fwd=config.defaults$seq_analysis$primer_action_fwd,
  primer_action_rev=config.defaults$seq_analysis$primer_action_rev,
  max_mismatches_fwd=config.defaults$seq_analysis$max_mismatches_fwd,
  max_mismatches_rev=config.defaults$seq_analysis$max_mismatches_rev,
  use_reverse_primers = config.defaults$seq_analysis$
    use_reverse_primers,
  reverse_primer_r1 = config.defaults$seq_analysis$reverse_primer_r1) {

  if (is.null(max_mismatches_fwd)) {
    max_mismatches_fwd <- max_mismatches
  }
  if (is.null(max_mismatches_rev)) {
    max_mismatches_rev <- max_mismatches
  }
  if (is.null(primer_action_fwd)) {
    primer_action_fwd <- primer_action
  }
  if (is.null(primer_action_rev)) {
    primer_action_rev <- primer_action
  }
  # Find matching loci, by primer sequence, for each of a set of sequences,
  # within a maximum distance.  The output rows will match the input sequences.
  match_primer_set <- function(
    seqs, locus_attrs, max_mismatches, colnm, prefix, do_revcmp) {
    primers <- locus_attrs[[colnm]]
    if (do_revcmp) {
      primers <- revcmp(primers)
    }
    matches <- find_primer_matches(seqs, primers, max_mismatches)
    # To ensure output and input match up, interpolate NA rows where needed
    matches <- matches[with(matches, order(SeqIdx, Mismatches)), ]
    matches <- matches[match(seq_along(seqs), matches$SeqIdx), ]
    matches$SeqIdx <- seq_along(seqs)
    # Assign a locus based on the primer match
    matches$Locus <- locus_attrs$Locus[matches$PrimerIdx]
    # Get rid of this index, since it only makes sense relative to the primer
    # data frame in here
    matches$PrimerIdx <- NULL
    # Add a custom prefix to each column
    mod <- colnames(matches) != "SeqIdx"
    colnames(matches)[mod] <- paste0(prefix, colnames(matches)[mod])
    matches
  }

  # Find a matching locus via forward primers
  matches_fwd <- match_primer_set(
    seqs, locus_attrs, max_mismatches, "Primer", "Fwd", FALSE)
  result <- matches_fwd
  # Same for reverse, if we're using those
  if (use_reverse_primers) {
    matches_rev <- match_primer_set(
      seqs, locus_attrs, max_mismatches, "ReversePrimer", "Rev",
      ! reverse_primer_r1)
    result <- merge(result, matches_rev, by.x = "SeqIdx", sort = FALSE)
    result$MatchingLocus <- result$FwdLocus
    result$MatchingLocus[with(
      result, is.na(FwdLocus) | is.na(RevLocus) | FwdLocus != RevLocus)] <- NA
  } else {
    result$MatchingLocus <- result$FwdLocus
  }
  result <- result[order(result$SeqIdx), ]
  result <- cbind(SeqOrig = seqs, result, stringsAsFactors = FALSE)
  result <- subset(result, select = -SeqIdx)
  result <- handle_primers(
    result, locus_attrs,
    primer_action_fwd, primer_action_rev, reverse_primer_r1)
  result
}

#' Modify reads based on matched primers
#' 
#' `handle_primers` carries out the read modification specified in
#' [make_read_primer_table].
#' 
#' @param result data frame of read and primer information as from
#'   [make_read_primer_table]
#' @param locus_attrs data frame of locus attributes
#' @param primer_action_fwd keyword for handling forward primers (see
#'   [make_read_primer_table])
#' @param primer_action_rev keyword for handling reverse primers (see
#'   [make_read_primer_table])
#' @param reverse_primer_r1 is each reverse primer in `locus_attrs` given in its
#'   orientation on the forward read?
#' @returns modified data frame with `Seq` column containing modified sequences
#' @md
handle_primers <- function(
  result, locus_attrs, primer_action_fwd, primer_action_rev, reverse_primer_r1) {
  # Handling this as two steps for both forward and reverse parts:
  #  1) substring to appropriate region
  #  2) add any replacement strings, for the special case of "replace"
  if (primer_action_rev != "none" && ! all(c("RevStart", "RevStop") %in% colnames(result))) {
    stop(
      paste0(
        "can't apply primer_action_rev \"",
        primer_action_rev,
        "\" without reverse primer information"))
  }
  start_poses <- switch(
    primer_action_fwd,
    none = 1,
    keep = result$FwdStart,
    replace = result$FwdStop + 1,
    remove = result$FwdStop + 1)
  if (is.null(start_poses)) {
    stop("primer_action_fwd should be one of none/keep/replace/remove")
  }
  end_poses <- switch(
    primer_action_rev,
    none = nchar(result$SeqOrig),
    keep = result$RevStop,
    replace = result$RevStart - 1,
    remove = result$RevStart - 1)
  if (is.null(end_poses)) {
    stop("primer_action_rev should be one of none/keep/replace/remove")
  }
  result$Seq <- substring(result$SeqOrig, start_poses, end_poses)
  if (primer_action_fwd == "replace") {
    prefix <- locus_attrs$Primer[
      match(result$MatchingLocus, locus_attrs$Locus)]
    result$Seq <- paste0(prefix, result$Seq)
  }
  if (primer_action_rev == "replace") {
    suffix <- locus_attrs$ReversePrimer[
      match(result$MatchingLocus, locus_attrs$Locus)]
    if (! reverse_primer_r1) {
      suffix <- revcmp(suffix)
    }
    result$Seq <- paste0(result$Seq, suffix)
  }
  idx <- is.na(result$MatchingLocus)
  result$Seq[idx] <- result$SeqOrig[idx]
  result
}

#' Find primer matches for reads
#' 
#' Returns a data frame of index values for the two input vectors for each
#' match, along with start and stop positions within each read and a count of
#' base mismatches.  There may be zero, one, or multiple output rows per input
#' read.
#' 
#' Indels are not supported.  Partial matches are supported at the 3' end and
#' are counted as mismatches.
#' 
#' If `max_mismatches` is `NA`, the best match for every read and primer
#' combination will be included in the output.  If `max_mismatches` is an
#' integer, at most one row will be provided for each read (for the first
#' discovered primer at or below that number of mismatches).  In this mode the
#' primers are checked for perfect matches first and then re-checked in
#' decreasing order of abundance to find any imperfect matches.  This means that
#' a search with a maximum mismatch count specified will be much faster than one
#' without, but could give counterintuitive results for some edge cases.  For
#' example, with max mismatches of 10 and a read containing primer A at 5
#' mismatches and primer B at 8 mismatches, either A or B may be matched
#' depending on the relative abundance of perfect matches found elsewhere. In
#' practice, with dissimilar primer sequences and a low number of allowed
#' mismatches, this is unlikely to be a problem.
#'
#' @param seqs_reads character vector of read sequences
#' @param seqs_primers character vector of primer sequences
#' @param max_mismatches integer number of mismatches allowed when checking
#'   primers against reads, or `NA` to check all combinations
#' @returns data frame of read and primer index pairs and match details
#' @md
find_primer_matches <- function(seqs_reads, seqs_primers, max_mismatches=NA) {
  # Any NAs will be handled in the same way as empty strings
  seqs_reads[is.na(seqs_reads)] <- ""
  seqs_primers[is.na(seqs_primers)] <- ""
  matches <- ordered(
    character(length(seqs_reads)),
    levels = unique(seqs_primers[seqs_primers != ""]))

  result_stub <- data.frame(
    SeqIdx = integer(),
    PrimerIdx = integer(),
    Start = integer(),
    Stop = integer(),
    Mismatches = integer())
  if (length(seqs_reads) == 0 || length(levels(matches)) == 0) {
    return(result_stub)
  }

  # Assign exact matches wherever observed, each time ignoring reads already
  # matched.  Using a munged version of the reads with some padding so we can
  # find partial matches around the edges (without dealing with actual gaps)
  if (! is.na(max_mismatches)) {
    padding <- paste0(rep("X", max(0, nchar(levels(matches)))), collapse = "")
    seqs_for_grep <- paste0(padding, seqs_reads, padding)
    for (seq_primer in levels(matches)) {
      idxl_mask <- is.na(matches)
      idx <- grep(seq_primer, seqs_for_grep[idxl_mask])
      matches[idxl_mask][idx] <- seq_primer
    }
    # reorder levels to prioritize the most abundant for all that follows
    matches <- stats::reorder(matches, -table(matches)[matches])
    # Second pass: approximate matches, in addition to those exact matches
    for (seq_primer in levels(matches)) {
      idxl_mask <- is.na(matches)
      idx <- agrep(
        seq_primer,
        seqs_for_grep[idxl_mask],
        max.distance = list(ins = 0, del = 0, all = max_mismatches))
      matches[idxl_mask][idx] <- seq_primer
    }
    matches <- stats::reorder(matches, -table(matches)[matches])
  }

  # Now check actual mismatch details
  len_max <- max(nchar(seqs_reads))
  # In this big matrix, each row is a position, each column a query; in the
  # colSums below, primers are compared as vectors column-wise to this matrix.
  # (Using raw bytes as that makes the comparisons significantly faster)
  queries <- make_raw_nt(seqs_reads, map = RAW_NT[1:4])
  slices <- seq_len(len_max)
  result <- do.call(rbind, lapply(levels(matches), function(seq_primer) {
    # Each per-primer data frame here will initially have a row for every read.
    # If max_mismatches is an integer we'll only assign to rows for reads
    # already determined to have a match for each primer, and then remove NA
    # rows in the end.  If max_mismatches is NA, we'll compare all to all and
    # keep all entries.
    # Empty primer sequences always yield NA rows out.
    result <- result_stub[seq_along(seqs_reads), ]
    result$SeqIdx <- seq_along(seqs_reads)
    result$PrimerIdx <- match(seq_primer, seqs_primers)
    primer <- make_raw_nt(seq_primer, pad = 0x40)[, 1]
    if (length(primer) > 1) {
      # If a max mismatch was given, we'll only check the reads that showed a
      # match for this primer.  Otherwise we'll check everything.
      if (is.na(max_mismatches)) {
        idxl_reads <- seq_along(seqs_reads)
        queries_here <- queries
      } else {
        idxl_reads <- matches %in% seq_primer
        queries_here <- queries[, idxl_reads, drop = FALSE]
      }
      if (ncol(queries_here) > 0) {
        for (start_pos in slices) {
          end_pos <- min(len_max, start_pos + nchar(seq_primer) - 1L)
          # to handle cases where the primer runs past the end of the reads,
          # we'll count those bases as mismatches
          truncate <- length(primer) - (end_pos - start_pos + 1)
          primer_segment <- primer[1:(end_pos - start_pos + 1)]
          chunk <- queries_here[start_pos:end_pos, , drop = FALSE]
          # If a pair of bytes have no overlapping bits, it's a mismatch.  Note
          # that also includes the 0x00 we're using as the padding value.
          mismatches <- truncate + colSums(
            (chunk & primer_segment) == as.raw(0))
          idxl_match <- with(
            result,
            is.na(Mismatches[idxl_reads]) | mismatches < Mismatches[idxl_reads])
          result$Start[idxl_reads][idxl_match] <- start_pos
          result$Stop[idxl_reads][idxl_match] <- end_pos
          result$Mismatches[idxl_reads][idxl_match] <- as.integer(mismatches[idxl_match])
        }
      }
    }
    subset(result, ! is.na(Mismatches))
  }))
  rownames(result) <- NULL
  result
}


# Util --------------------------------------------------------------------

#' IUPAC DNA nucleotides as raw bytes
#' 
#' `RAW_NT` is a named vector of raw bytes allowing fast comparison of IUPAC
#' DNA nucleotide characters.  A, C, G, and T are assigned single individual
#' bits, and the standard IUPAC characters are assigned bitwise combinations of
#' these.  This way `RAW_NT["A"] & RAW_NT["G"]` is zero, while
#' `RAW_NT["A"] & RAW_NT["R"]` is nonzero.  The gap characters - and . are
#' assigned one remaining bit, leaving 3 bits unspecified.
#' 
#' ```
#'      - TGCA
#' A 0000 0001  01
#' C 0000 0010  02
#' G 0000 0100  04
#' T 0000 1000  08
#' 
#' R 0000 0101  05
#' Y 0000 1010  0A
#' S 0000 0110  06
#' W 0000 1001  09
#' K 0000 1100  0C
#' M 0000 0011  03
#' 
#' B 0000 1110  0E
#' D 0000 1101  0D
#' H 0000 1011  0B
#' V 0000 0111  07
#' 
#' N 0000 1111  0F
#' 
#' - 0001 0000  10
#' . 0001 0000  10
#' ```
#' @md
RAW_NT <- as.raw(c(1, 2, 4, 8, 16, 16))
names(RAW_NT) <- c("A", "C", "G", "T", "-", ".")
RAW_NT["R"] <- RAW_NT["A"] | RAW_NT["G"]
RAW_NT["Y"] <- RAW_NT["C"] | RAW_NT["T"]
RAW_NT["S"] <- RAW_NT["G"] | RAW_NT["C"]
RAW_NT["W"] <- RAW_NT["A"] | RAW_NT["T"]
RAW_NT["K"] <- RAW_NT["G"] | RAW_NT["T"]
RAW_NT["M"] <- RAW_NT["A"] | RAW_NT["C"]
RAW_NT["B"] <- RAW_NT["C"] | RAW_NT["G"] | RAW_NT["T"]
RAW_NT["D"] <- RAW_NT["A"] | RAW_NT["G"] | RAW_NT["T"]
RAW_NT["H"] <- RAW_NT["A"] | RAW_NT["C"] | RAW_NT["T"]
RAW_NT["V"] <- RAW_NT["A"] | RAW_NT["C"] | RAW_NT["G"]
RAW_NT["N"] <- RAW_NT["A"] | RAW_NT["C"] | RAW_NT["G"] | RAW_NT["T"]

#' Complements of IUPAC nucleotide codes
#' 
#' `CMP` is a named character vector mapping each IUPAC DNA nucleotide code to
#' its complement.  For example, `CMP["A"]` is `T`.
#' @md
CMP <- c(
  A = "T",
  C = "G",
  T = "A",
  G = "C",
  R = "Y",
  Y = "R",
  S = "S",
  W = "W",
  K = "M",
  M = "K",
  B = "V",
  D = "H",
  H = "D",
  V = "B")

#' Make matrix of raw bytes from nucleotide sequences
#' 
#' Each sequence in the given vector becomes a column of the output matrix, with
#' a row for each position.  Shorter sequences are padded with a specific value
#' at bottom of the matrix.  With the defaults (see [RAW_NT]), A, C, G, and T
#' (case insensitive) are encoded as 01, 02, 03, and 04, IUPAC codes are bitwise
#' combinations of those values, padding values are 0x80, and any other
#' character is 0x00.
#' 
#' @param seqs character vector of nucleotide sequences
#' @param map raw vector with nucleotide names and byte values
#' @param pad raw value to use for missing positions
#' @param other raw value to use for nucleotides not in `map`
#' @return raw matrix with positions on rows and sequences on columns
#' @md
make_raw_nt <- function(seqs, map = RAW_NT, pad = 0x80, other = 0x00) {
  pad <- as.raw(as.integer(pad))
  other <- as.raw(as.integer(other))
  seq_levels <- names(map)
  if (length(seqs) == 0) {
    return(matrix(raw())[0, 0])
  }
  seqs <- toupper(seqs)
  len <- max(nchar(seqs))
  out <- do.call(cbind, lapply(strsplit(seqs, ""), function(vec){
    # will put NA for too-short seqs, but we need to distinguish that from NA
    # from unrecognized NTs
    vec <- vec[seq_len(len)]
    pads <- is.na(vec)
    # First set all NA entries to "other" value, but mask the missing cases with
    # the pad value
    vec <- map[vec]
    vec[is.na(names(vec))] <- other
    vec[pads] <- pad
    unname(vec)
  }))
  out
}
