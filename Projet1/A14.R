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


############
### A1.3 ###
############

dimers_frequency <- function(sequence){
  #'Function to compute de frequency of each dimer along a sequence
  #'
  #'    @param sequence : a vector of characters corresponding to a nucleotide sequence (sequence of capital "A"s, "T"s, "C"s and, "G"s)
  #'    
  #'    @return A list containing: 
  #'                An RDS matrix containing the frequency of each dimer and save it as "observed_frequency.rds"
  #'                A vector containing the under represented dimers, ordered increasingly
  #'                A vector containing the over represented dimers, ordered increasingly
  #'
  dimers <- sliding_windows(sequence = sequence, size = 2, shift = 1)
  
  col_names <- c("*A", "*C", "*G", "*T")
  row_names <- c("A*", "C*", "G*", "T*")
  
  observed_dimers <- matrix(nrow = 4, ncol = 4)
  rownames(observed_dimers) <- row_names
  colnames(observed_dimers) <- col_names
  
  under_represented <- list()
  over_represented <- list()
  
  for (i in seq_along(row_names)){
    for (j in seq_along(col_names)){
      dimer <- paste0(strsplit(row_names[[i]], split = "")[[1]][1], strsplit(col_names[[j]], split = "")[[1]][2])
      count <- sum(dimers == dimer)
      observed_dimers[i, j] <- round(count / length(dimers), digits=4)
      if (round(count / length(dimers), digits=4) < 1/16) {
        under_represented[[length(under_represented) + 1]] <- c(dimer, round(count / length(dimers), digits=4))
      } 
      if (round(count / length(dimers), digits=4) > 1/16) {
        over_represented[[length(over_represented) + 1]] <- c(dimer, round(count / length(dimers), digits=4))
      }
    }
  }
  
  under_repr_order <- order(sapply(under_represented, function(x) as.numeric(x[2])))
  sorted_under_repr <- under_represented[under_repr_order]
  over_repr_order <- order(sapply(over_represented, function(x) as.numeric(x[2])))
  sorted_over_repr <- over_represented[over_repr_order]
  
  saveRDS(observed_dimers, file = "observed_frequency.rds")
  return(list(observed_dimers=observed_dimers, under_represented=sorted_under_repr, over_represented=sorted_over_repr))
}


################################ End of previous functions ################################ 




cat("Getting the Bovine papular stomatitis virus ... \n")
sequence <- toupper(read.GenBank(access.nb = "NC_005337", as.character = TRUE)[[1]])
cat("The lenght of the Bovine papular stomatitis virus is: ", length(sequence))


odd_ratios <- function(sequence, observed_frequencies, k){
  #'Function to compute the odd ratio of k-mers from a sequence
  #'
  #'    @param sequence : A vector of characters corresponding to a nucleotide sequence (sequence of capital "A"s, "T"s, "C"s and, "G"s)
  #'    @param table : An RDS matrix where each cell is the observed frequency for the k-mer
  #'    @param k : A numerical value that indicate the k-mer (k=2 for dimer, 3 for 3-mer, ...)
  #'    
  #'    @return A list containing: 
  #'        An RDS matrix of the same size as table with the odd ratio.
  #'        A vector containing the under represented dimers, ordered increasingly
  #'        A vector containing the over represented dimers, ordered increasingly
  
  expected_frequency <- c(
                          A = (sum(sequence == "A")/length(sequence)), 
                          T = (sum(sequence == "T")/length(sequence)), 
                          G = (sum(sequence == "G")/length(sequence)), 
                          C = (sum(sequence == "C")/length(sequence))
                        )
  odd_ratio_table <- observed_frequencies
  
  represented <- list()

  for (i in seq_along(rownames(odd_ratio_table))){
    for (j in seq_along(colnames(odd_ratio_table))){
      row_ntd <- strsplit(rownames(odd_ratio_table)[[i]], split = "")[[1]][1]
      col_ntd <- strsplit(colnames(odd_ratio_table)[[j]], split = "")[[1]][2]
      dimer <- paste0(row_ntd, col_ntd)
      odd_ratio_table[i, j] <- round(observed_frequencies[i, j] / (expected_frequency[row_ntd] * expected_frequency[col_ntd]), digits = 4)
      
      represented[[length(represented) + 1]] <- c(dimer, round(odd_ratio_table[i, j], digits=4))
    }
  }
  
  repr_order <- order(sapply(represented, function(x) as.numeric(x[2])))
  sorted_repr <- represented[repr_order]
  saveRDS(odd_ratio_table, file = "odd_ratio_table.rds")
  return(list(odd_ratio_table=odd_ratio_table, sorted_repr=sorted_repr))
}

frequency <- dimers_frequency(sequence)
ratios <- odd_ratios(sequence, frequency$observed_dimers, 2)










