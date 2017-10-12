# simulated data for testing

txt.locus_attrs <- textConnection(
"Locus   LengthMin  LengthMax    LengthBuffer   Motif   Primer                  ReversePrimer
A        131        179          20             TAGA    TATCACTGGTGTTAGTCCTCTG  CACAGTTGTGTGAGCCAGTC
B        194        235          20             TAGA    AGTCTCTCTTTCTCCTTGCA    TAGGAGCCTGTGGTCCTGTT
1        232        270          20             TATC    ACAGTCAAGAATAACTGCCC    CTGTGGCTCAAAAGCTGAAT
2        218        337          20             TCCA    TTGTCTCCCCAGTTGCTA      TCTGTCATAAACCGTCTGCA
")
locus_attrs <- read.table(txt.locus_attrs, header = T, row.names = 1, stringsAsFactors = F)
close(txt.locus_attrs)

make.seq_junk <- function(N) {
  nucleotides <- c("A", "T", "C", "G")
  vapply(runif(N, min = 1, max = 20), function(L)
    paste0(sample(nucleotides, L, replace = T), collapse = ""),
    "character")
}

simulate.seqs <- function(locus_name, locus_attrs, homozygous=NULL, N=5000,
                          off_target_ratio=10) {
  attrs <- locus_attrs[locus_name, ]
  L.min <- attrs$LengthMin - nchar(attrs$Primer)
  L.max <- attrs$LengthMax - nchar(attrs$Primer)
  L <- as.integer(runif(2, min = L.min, max = L.max))
  if (!missing(homozygous)) {
    if (homozygous) {
      L[2] <- L[1]
    } else {
      while (L[1] == L[2])
      L <- as.integer(runif(2, min = L.min, max = L.max))
    }
  }
  make.reps <- function(motif, len)
    paste(rep(motif, floor(len / nchar(motif))), collapse = "")
  repeats1 <- make.reps(attrs$Motif, L[1])
  repeats2 <- make.reps(attrs$Motif, L[2])
  seq.correct.1 <- paste0(attrs$Primer, repeats1, attrs$ReversePrimer)
  seq.correct.2 <- paste0(attrs$Primer, repeats2, attrs$ReversePrimer)
  # The exact correct sequences, either one or two unique
  N1 <- N / 2 + rnorm(10000) * N / 20
  N1 <- as.integer(max(min(N1, N * 0.9), N * 0.1))
  N2 <- N - N1
  seqs <- c(rep(seq.correct.1, N1),
            rep(seq.correct.2, N2))
  ## Mix in stutter
  N1.stutter <- as.integer(runif(1, min = N1 / 20, max = N1 / 3))
  N2.stutter <- as.integer(runif(1, min = N2 / 20, max = N2 / 3))
  seqs[1:N1.stutter] <- paste0(attrs$Primer,
                               make.reps(attrs$Motif,
                                         L[1] - nchar(attrs$Motif)),
                               attrs$ReversePrimer)
  seqs[N1 + 1:N2.stutter] <- paste0(attrs$Primer,
                                  make.reps(attrs$Motif,
                                            L[2] - nchar(attrs$Motif)),
                                  attrs$ReversePrimer)
  ## Mix in off-target amplification
  idx <- seq(1, length(seqs), off_target_ratio)
  seqs[idx] <- paste0(attrs$Primer, make.seq_junk(10), attrs$ReversePrimer)
  ## TODO add contam
  ## TODO munge up the reads a bit
  ## TODO shuffle vector
  return(seqs)
}

simulate.set <- function(locus_attrs) {
  seqs <- lapply(rownames(locus_attrs), function(n) {
    s <- simulate.seqs(n, locus_attrs)
    s
  })
  names(seqs) <- rownames(locus_attrs)
  return(seqs)
}

# setting the set to an arbitrary value so the below is all arbitrary but
# still deterministic.
set.seed(0)
# Three sets of samples across all loci
seqs1 <- simulate.set(locus_attrs)
seqs2 <- simulate.set(locus_attrs)
seqs3 <- simulate.set(locus_attrs)
set.seed(NULL)
