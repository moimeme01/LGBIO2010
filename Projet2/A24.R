library(Biostrings)
library(pwalign)
library(ape)

# ==============================================================================
# 1. Function to compute the relative substitution (mismatch) frequency matrix
# ==============================================================================
computeRelativeFrequencies <- function(seqs, substitutionMatrix, gapOpening, gapExtension) {
  n <- length(seqs)
  seq_names <- names(seqs)
  seq_chars <- as.character(seqs)

  # Initialize n x n matrix with zeros (guarantees m[i, i] = 0)
  m <- matrix(0, nrow = n, ncol = n, dimnames = list(seq_names, seq_names))

  # Compute only the upper triangle (i < j) and mirror it so m[i, j] == m[j, i] strictly
  if (n > 1) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        align <- pairwiseAlignment(
          pattern = seq_chars[i],
          subject = seq_chars[j],
          type = "global",
          substitutionMatrix = substitutionMatrix,
          gapOpening = gapOpening,
          gapExtension = gapExtension
        )

        n_match <- as.numeric(nmatch(align))
        n_mismatch <- as.numeric(nmismatch(align))
        denom <- n_match + n_mismatch

        rel_freq <- if (denom > 0) n_mismatch / denom else 0
        m[i, j] <- rel_freq
        m[j, i] <- rel_freq
      }
    }
  }

  return(m)
}

# ==============================================================================
# 2. Build StudySet: [2.1] SARS-CoV-2 subset + 3 sequences from [2.3] Q2
# ==============================================================================
# Load the subset selected in task [2.1]
sars2_seqs <- readDNAStringSet("./Projet2/FilteredSequences.fasta")

# Ensure names are clean AC numbers (strip any FASTA description after a space)
names(sars2_seqs) <- sub(" .*", "", names(sars2_seqs))

# Fetch the 3 potentially related sequences from [2.3] Question 2
related_acs <- c("Y10051.1", "CY125943.1", "A22884.1")
raw_related <- read.GenBank(related_acs, species.names = FALSE, as.character = TRUE)
related_strings <- vapply(
  raw_related,
  function(x) paste0(toupper(x), collapse = ""),
  FUN.VALUE = character(1)
)
related_seqs <- DNAStringSet(related_strings)
names(related_seqs) <- related_acs

# Combine into StudySet (avoiding duplicates if re-run)
StudySet <- c(sars2_seqs, related_seqs[!names(related_seqs) %in% names(sars2_seqs)])

# Optional: save StudySet locally so you don't have to re-download it in later questions
saveRDS(StudySet, file = "StudySet.rds")

# ==============================================================================
# 3. Compute the matrix using Question 4 parameters from task [2.2]
# ==============================================================================
substitutionMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
gapOpening <- 5
gapExtension <- 2

m <- computeRelativeFrequencies(
  seqs = StudySet,
  substitutionMatrix = substitutionMatrix,
  gapOpening = gapOpening,
  gapExtension = gapExtension
)

# ==============================================================================
# 4. Verify matrix properties and save to .rds for upload
# ==============================================================================
stopifnot(
  is.matrix(m),
  all(diag(m) == 0),
  all(m >= 0),
  isSymmetric(m),
  identical(rownames(m), colnames(m)),
  !any(is.na(rownames(m)))
)

print(round(m, 4))
cat("Sequence labels (AC numbers):\n")
print(rownames(m))

# Save the matrix variable to an .rds file in your working directory
saveRDS(m, file = "relative_frequencies.rds")
cat("Saved matrix successfully to:", normalizePath("relative_frequencies.rds"), "\n")

#q2 
freq2dist = function(freq, nbStates) {
  ratio <- (nbStates - 1) / nbStates
  dist <- -ratio * log(1 - freq / ratio)
  return(dist + 0)
}

#q3

freq2dist = function(freq, nbStates) {
  ratio <- (nbStates - 1) / nbStates
  return(-ratio * log(1 - freq / ratio) + 0)
}

m <- readRDS("./relative_frequencies.rds")
dist_mat <- freq2dist(m, nbStates = 4)

cat("Maximal JC evolutionary distance:", max(dist_mat), "\n")