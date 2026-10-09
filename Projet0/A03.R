library(ape)

seq_length <- function(fasta_file) {
  sequence <- read.dna(fasta_file, format = "fasta", as.character = TRUE)
  len <- length(sequence)
  len
}

seq_length("NC_002662.fasta")

