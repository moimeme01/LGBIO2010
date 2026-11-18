# A2.1

library(ape)
#install.packages("seqinr")
library(seqinr)
library(Biostrings)

# Accession numbers
accession_numbers <- scan("/Users/thibaultvanni/Documents/UCL/LGBIO2010 - Bioinformatics/LGBIO2010-Projet2/AccessionNumbers.txt", what = character(), quiet = TRUE)

accession_numbers <- sort(unique(accession_numbers))

# Get sequences
seq_list <- read.GenBank(accession_numbers, species.names = FALSE, as.character = TRUE) 

##################
### Question 1 ###
##################

#1# Keep only the found ones and puts in increasing AC order
if (length(seq_list) == 0) {
  cat(0, "\n")
  quit(save = "no")
}

seq_list <- seq_list[sort(names(seq_list))]


# Convert to uppercase strings 
seq_strings <- vapply(
  seq_list,
  function(x) toupper(c2s(getSequence(x))),
  FUN.VALUE = character(1)
)


#2# Keep sequences made exclusively of A/T/C/G 
is_atcg <- grepl("^[ATCG]+$", seq_strings)
seq_strings <- seq_strings[is_atcg]
#print(seq_strings)

# If nothing remains after filtering, print 0 and exit
if (length(seq_strings) == 0) {
  cat(0, "\n")
  quit(save = "no")
}



#3# Remove duplicates by exact nucleotide content keeping the first 
keep <- !duplicated(seq_strings)
seq_strings_unique <- seq_strings[keep]
length <- length(seq_strings_unique)

seq_set <- DNAStringSet(seq_strings_unique)
writeXStringSet(seq_set, filepath = "/Users/thibaultvanni/Documents/UCL/LGBIO2010 - Bioinformatics/LGBIO2010-Projet2/FilteredSequences.fasta")

print("The number of sequences is :")
print(length)

##################
### Question 2 ###
##################

avg_length <- mean(nchar(seq_strings_unique)) 
print("The average sequence length")
print(avg_length)  

