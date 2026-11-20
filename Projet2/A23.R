library(Biostrings)
library(pwalign)
library(ape)

pairwiseAlignmentResult = function(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension) {
  align <- pairwiseAlignment(
    pattern = pattern,
    subject = subject,
    substitutionMatrix = substitutionMatrix,
    gapOpening = gapOpening,
    gapExtension = gapExtension,
    type = type
  )
  totalScore <- score(align)
  numbmatch <- as.numeric(nmatch(align))
  numbmismatch <- as.numeric(nmismatch(align))
  indelObj <- nindel(align)
  nGap <- sum(indelObj@insertion[, "WidthSum"]) + sum(indelObj@deletion[, "WidthSum"])

  alignmentLength <- numbmatch + numbmismatch + nGap
  normalisedScore <- totalScore / alignmentLength
  return(normalisedScore)
}

# q1
pairwiseSignificance = function(alpha, pattern, subject, type, substitutionMatrix, gapOpening, gapExtension) {
  getScore <- function(res) {
    if (is.numeric(res)) {
      return(as.numeric(res)[1])
    } else if (is.list(res)) {
      if (!is.null(res$score)) return(as.numeric(res$score)[1])
      num_elems <- Filter(is.numeric, res)
      if (length(num_elems) > 0) return(as.numeric(num_elems[[1]])[1])
    } else if (isS4(res)) {
      if (exists("score", mode = "function")) return(as.numeric(score(res))[1])
      if ("score" %in% slotNames(res)) return(as.numeric(res@score)[1])
    }
    stop("Unable to extract alignment score from pairwiseAlignmentResult().")
  }

  obs_score <- getScore(
    pairwiseAlignmentResult(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension)
  )

  n <- max(1000L, as.integer(ceiling(10 / alpha)))
  pattern_chars <- strsplit(pattern, "", fixed = TRUE)[[1]]
  subject_chars <- strsplit(subject, "", fixed = TRUE)[[1]]

  random_scores <- numeric(n)
  for (i in seq_len(n)) {
    rand_pattern <- paste0(sample(pattern_chars), collapse = "")
    rand_subject <- paste0(sample(subject_chars), collapse = "")
    random_scores[i] <- getScore(
      pairwiseAlignmentResult(rand_pattern, rand_subject, type, substitutionMatrix, gapOpening, gapExtension)
    )
  }

  threshold <- quantile(random_scores, probs = 1 - alpha, names = FALSE)
  return(as.logical(obs_score > threshold))
}

# q2
substitutionMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
gapOpening <- 5
gapExtension <- 2

accessions <- c("MW090859.1", "Y10051.1", "CY125943.1", "A22884.1")
raw_seqs <- read.GenBank(accessions, species.names = FALSE, as.character = TRUE)
dna_seqs <- vapply(
  raw_seqs,
  function(x) paste0(toupper(x), collapse = ""),
  FUN.VALUE = character(1)
)

set.seed(42)
sig_results <- vapply(
  dna_seqs[c("Y10051.1", "CY125943.1", "A22884.1")],
  function(s) {
    pairwiseSignificance(
      alpha = 0.05,
      pattern = dna_seqs[["MW090859.1"]],
      subject = s,
      type = "global",
      substitutionMatrix = substitutionMatrix,
      gapOpening = gapOpening,
      gapExtension = gapExtension
    )
  },
  FUN.VALUE = logical(1)
)

print(sig_results)
cat("number of significant alignments:", sum(sig_results), "\n")