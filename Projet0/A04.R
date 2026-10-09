library(ape)

nb_sliding_windows <- function(l, size, shift) {
  count <- 0
  position <- 0
  while (position + size < l) {
    count <- count + 1
    position <- position + shift
  }
  count <- count + 1
  count
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


frequency_in_sliding_window <- function(sequence, ntd, size, shift) {
  nb <- nb_sliding_windows(length(sequence), size, shift)
  up_sequence <- toupper(sequence)
  
  resultat <- numeric(nb) # We need to do a vector because we will replace each value by a list

  windows <- sliding_windows(sequence, size, shift)
  
  index <- 1
  for (element in windows){
    cat("element is: ", element, "\n")
    count <- 0
    for (letter in strsplit(element, split='')[[1]]){
      if (letter == ntd) {
        count <- count + 1
      }
    }
    resultat[[index]] = count/length(strsplit(element, split='')[[1]])
    index <- index + 1
  }
  
  return (resultat)
}


seq <- c("C", "A", "T", "A", "C", "C", "T", "T", "C", "C", "A", "A", "A")

frequency_in_sliding_window(seq, "T", 6, 6)


