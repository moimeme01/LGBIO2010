# Make sure Biostrings and pwalign are both installed and loaded
library(Biostrings)
library(pwalign)
library(ape)

######################
### Question 1 - 2 ###
######################

# Get a valid DNA substitution matrix (for match = +2, mismatch = -1)
submat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

# Perform the global alignment with linear gap penalty (d = 4)
align <- pairwiseAlignment(
  pattern = "ATTCCGA",
  subject = "TATCG",
  substitutionMatrix = submat,
  gapOpening = 0,
  gapExtension = 4,
  type = "overlap"        #Changer en global pour la question 1 et overlap pour la question 2
)

totalScore <- score(align)
numbmatch <- as.numeric(nmatch(align))
numbmismatch <- as.numeric(nmismatch(align))
# nindel may return InDel object
indelObj <- nindel(align)

# Sum width of insertions and deletions
nGap <- sum(indelObj@insertion[,"WidthSum"]) + sum(indelObj@deletion[,"WidthSum"])



##################
### Question 3 ###
##################

pairwiseAlignmentResult = function(pattern, subject, type, substitutionMatrix, gapOpening, gapExtension) {
  align <- pairwiseAlignment(pattern = pattern, subject = subject, substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension, type = type)
  totalScore <- score(align)
  numbmatch <- as.numeric(nmatch(align))
  numbmismatch <- as.numeric(nmismatch(align))
  # nindel may return InDel object
  indelObj <- nindel(align)
  
  # Sum width of insertions and deletions
  nGap <- sum(indelObj@insertion[,"WidthSum"]) + sum(indelObj@deletion[,"WidthSum"])
  
  
  alignmentLength <- numbmatch + numbmismatch + nGap
  
  normalisedScore <- totalScore / alignmentLength
  
  
  return(normalisedScore)
}



test <- pairwiseAlignmentResult(pattern = "ATTCCGA", subject = "ATCGA", type = "global", substitutionMatrix = submat, gapOpening = 0, gapExtension = 4)

print(test)




##################
### Question 4 ###
##################

substitutionMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
gapOpening <- 5
gapExtension <- 2


fasta_file <- "./FilteredSequences.fasta"
seqs <- readDNAStringSet(fasta_file)

toCompare <- read.GenBank("MW090859.1", species.names = FALSE, as.character = TRUE) 
toCompare <- vapply(
  toCompare,
  function(x) toupper(c2s(getSequence(x))),
  FUN.VALUE = character(1)
)


scores <- vapply(
  seqs,
  function(s) pairwiseAlignmentResult(as.character(s), as.character(toCompare), type = "global", substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension),
  FUN.VALUE = numeric(1)
)

normalised <- sum(scores)/length(scores)
print("The average alignment score is")
print(normalised)

##################
### Question 5 ###
##################



substitutionMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
gapOpening <- 10000
gapExtension <- 2


fasta_file <- "./FilteredSequences.fasta"
seqs <- readDNAStringSet(fasta_file)

toCompare <- read.GenBank("MW090859.1", species.names = FALSE, as.character = TRUE) 
toCompare <- vapply(
  toCompare,
  function(x) toupper(c2s(getSequence(x))),
  FUN.VALUE = character(1)
)


scores <- vapply(
  seqs,
  function(s) pairwiseAlignmentResult(as.character(s), as.character(toCompare), type = "global", substitutionMatrix = substitutionMatrix, gapOpening = gapOpening, gapExtension = gapExtension),
  FUN.VALUE = numeric(1)
)

normalised <- sum(scores)/length(scores)

print(normalised)



