library(ape)


####################################
### Previously created functions ###
####################################

############
### A0.4 ###
############

nb_sliding_windows <- function(l, size, shift) {
  #'Function to compute the number of sliding window positions
  #'    @param l : a numeric value representing the length, in base pairs, of a sequence,
  #'    @param size: a numeric value representing the size of the sliding window,
  #'    @param shift : a numeric value representing the shift between successive sliding window positions,
  #'
  #'    @return numeric value, which is the number of possible sliding window positions given the above parameters
  
  count <- 0
  position <- 0
  while (position + size < l) {
    count <- count + 1
    position <- position + shift
  }
  count <- count + 1
  return(count)
}

sliding_windows <- function(sequence, size, shift){
  #'Function to create the windows 
  #'    @param sequence : a vector containing the nucleotides
  #'    @param size : the size of the windows
  #'    @param shift : the number of nucloetides shifted between windows
  #'
  #'    @return a vector containing strings of nucleotides of length size 
  
  nb <- nb_sliding_windows(length(sequence), size, shift)
  resultat <- vector(mode="list", length = nb)
  str_seq <- toupper(sequence)
  
  start <- 1
  for (i in seq_along(resultat)){
    end <- min(start + size - 1, length(str_seq))
    resultat[[i]] <- paste0(str_seq[start:end], collapse = "")
    start <- start + shift
  }
  return (resultat)
}


############
### A1.1 ###
############


gc_content <- function(sequence, size, shift){
  #' Compute the GC-content of a sequence using sliding window.
  #' 
  #'    @param sequence: a vector of characters corresponding to a nucleotide sequence (sequence of capital "A"s, "T"s, "C"s and, "G"s)
  #'    @param size: a numeric value corresponding to the size of the sliding window
  #'    @param shift: a numeric value representing the shift between successive sliding window positions
  #'    
  #'@return a numeric vector where the i-th value corresponds to the GC-content of the sequence computed on the i-th position of the sliding window.
  nb_windows <- nb_sliding_windows(length(sequence), size, shift)
  
  resultat <- numeric(nb_windows)
  windows <- sliding_windows(sequence, size, shift)
  
  for (i in seq_along(resultat)){
    cat("considering: ", windows[[i]], "\n")
    letters <- strsplit(windows[[i]], split="")[[1]]
    count <- sum(letters == "G") + sum(letters == "C")
    ratio <- count/length(letters)
    resultat[i] <- ratio
  }
  return (resultat)
}


################################ End of previous functions ################################ 


sequence <- read.dna(file="sequence.fasta", format = "fasta")
cat("The length of the sequence is: ", length(sequence), "\n")

