library(ape)

# Téléchargement de la séquence depuis GenBank
seq_sp <- read.GenBank(access.nb = "NZ_JFIS01000023.1", as.character = TRUE)[[1]]

longueur_bp <- length(seq_sp)
cat("long seq :", longueur_bp, "pb\n")

library(ape)
sequence <- read.GenBank(access.nb = "NZ_JFIS01000023.1", as.character = TRUE)[[1]]

set.seed(42)
null_sequence <- sample(sequence)

bacterial_starts <- c("TTG", "CTG", "ATA", "ATT", "ATC", "ATG", "GTG")
standard_stops   <- c("TGA", "TAA", "TAG")
null_orfs_lengths <- find_all_orfs(null_sequence, bacterial_starts, standard_stops)

max_null_length <- max(null_orfs_lengths)
cat("long max d'un orf null", max_null_length, "pb\n")


real_orfs_lengths <- find_all_orfs(sequence, bacterial_starts, standard_stops)
nb_candidate_genes <- sum(real_orfs_lengths > max_null_length)

cat("orf plus grand que de base ", nb_candidate_genes, "\n")
threshold_p01 <- quantile(null_orfs_lengths, 0.99)
nb_candidates_p01 <- sum(real_orfs_lengths > threshold_p01)
cat("seuil longueur :", threshold_p01, "pb\n")
cat("nb gene candidats :", nb_candidates_p01, "\n")


find_longest_candidate_frame <- function(sequence, starts = bacterial_starts, stops = standard_stops) {
  seq_upper <- toupper(sequence)
  if (length(seq_upper) == 1) {
    seq_upper <- strsplit(seq_upper, "")[[1]]
  }
  L <- length(seq_upper)
  
  best_orf <- list(length = 0, frame = "", strand = "", start = 0, end = 0)
  
  for (f in 0:2) {
    sub_seq <- seq_upper[(1 + f):L]
    n_codons <- floor(length(sub_seq) / 3)
    if (n_codons == 0) next
    
    idx_start <- seq(1, n_codons * 3, by = 3)
    codons <- paste0(sub_seq[idx_start], sub_seq[idx_start + 1], sub_seq[idx_start + 2])
    
    in_orf <- FALSE
    start_pos <- 0
    
    for (i in seq_along(codons)) {
      codon <- codons[i]
      if (!in_orf) {
        if (codon %in% starts) {
          in_orf <- TRUE
          start_pos <- i
        }
      } else {
        if (codon %in% stops) {
          len_nt <- (i - start_pos) * 3
          if (len_nt > best_orf$length) {
            start_genomic <- 1 + f + (start_pos - 1) * 3
            frame_val <- paste0("+", (start_genomic - 1) %% 3 + 1)
            
            best_orf <- list(
              length = len_nt,
              frame = frame_val,
              strand = "+",
              start = start_genomic,
              end = 1 + f + i * 3 - 1
            )
          }
          in_orf <- FALSE
        }
      }
    }
  }
  
  rev_seq <- reverse_complement(seq_upper)
  for (f in 0:2) {
    sub_seq <- rev_seq[(1 + f):L]
    n_codons <- floor(length(sub_seq) / 3)
    if (n_codons == 0) next
    
    idx_start <- seq(1, n_codons * 3, by = 3)
    codons <- paste0(sub_seq[idx_start], sub_seq[idx_start + 1], sub_seq[idx_start + 2])
    
    in_orf <- FALSE
    start_pos <- 0
    
    for (i in seq_along(codons)) {
      codon <- codons[i]
      if (!in_orf) {
        if (codon %in% starts) {
          in_orf <- TRUE
          start_pos <- i
        }
      } else {
        if (codon %in% stops) {
          len_nt <- (i - start_pos) * 3
          if (len_nt > best_orf$length) {
            rev_start_pos <- 1 + f + (start_pos - 1) * 3
            orig_5prime <- L - rev_start_pos + 1
            frame_num <- - ((orig_5prime %% 3) + 1)
            
            best_orf <- list(
              length = len_nt,
              frame = as.character(frame_num),
              strand = "-",
              rev_start = rev_start_pos,
              orig_5prime = orig_5prime
            )
          }
          in_orf <- FALSE
        }
      }
    }
  }
  
  return(best_orf)
}

longest_orf_res <- find_longest_candidate_frame(sequence)

cat("long pb :", longest_orf_res$length, "\n")
cat("Brin:", longest_orf_res$strand, "\n")
cat("Cadre de lec", longest_orf_res$frame, "\n")



cat("pos de base ", longest_orf_res$start, "\n")




################# version clean pour ingi




library(ape)

bacterial_starts <- c("TTG", "CTG", "ATA", "ATT", "ATC", "ATG", "GTG")
standard_stops <- c("TGA", "TAA", "TAG")

reverse_complement <- function(seq_vector) {
  chartr("ATGCatgc", "TACGtacg", rev(seq_vector))
}

find_orfs_strand <- function(sequence, starts, stops) {
  orfs_lengths <- numeric()
  seq_upper <- toupper(sequence)
  seq_len <- length(seq_upper)
  
  for (frame in 0:2) {
    sub_seq <- seq_upper[(1 + frame):seq_len]
    n_codons <- floor(length(sub_seq) / 3)
    if (n_codons == 0) next
    
    idx_start <- seq(1, n_codons * 3, by = 3)
    codons <- paste0(sub_seq[idx_start], sub_seq[idx_start + 1], sub_seq[idx_start + 2])
    
    in_orf <- FALSE
    start_pos <- 0
    
    for (i in seq_along(codons)) {
      codon <- codons[i]
      if (!in_orf) {
        if (codon %in% starts) {
          in_orf <- TRUE
          start_pos <- i
        }
      } else {
        if (codon %in% stops) {
          len_nt <- (i - start_pos) * 3
          orfs_lengths <- c(orfs_lengths, len_nt)
          in_orf <- FALSE
          start_pos <- 0
        }
      }
    }
  }
  return(orfs_lengths)
}

find_all_orfs <- function(sequence, starts = bacterial_starts, stops = standard_stops) {
  if (length(sequence) == 1) {
    sequence <- strsplit(sequence, "")[[1]]
  }
  fwd_orfs <- find_orfs_strand(sequence, starts, stops)
  rev_seq <- reverse_complement(sequence)
  rev_orfs <- find_orfs_strand(rev_seq, starts, stops)
  return(c(fwd_orfs, rev_orfs))
}

find_longest_candidate_frame <- function(sequence, starts = bacterial_starts, stops = standard_stops) {
  seq_upper <- toupper(sequence)
  if (length(seq_upper) == 1) {
    seq_upper <- strsplit(seq_upper, "")[[1]]
  }
  L <- length(seq_upper)
  
  best_orf <- list(length = 0, frame = "", strand = "", start = 0, end = 0)
  
  for (f in 0:2) {
    sub_seq <- seq_upper[(1 + f):L]
    n_codons <- floor(length(sub_seq) / 3)
    if (n_codons == 0) next
    
    idx_start <- seq(1, n_codons * 3, by = 3)
    codons <- paste0(sub_seq[idx_start], sub_seq[idx_start + 1], sub_seq[idx_start + 2])
    
    in_orf <- FALSE
    start_pos <- 0
    
    for (i in seq_along(codons)) {
      codon <- codons[i]
      if (!in_orf) {
        if (codon %in% starts) {
          in_orf <- TRUE
          start_pos <- i
        }
      } else {
        if (codon %in% stops) {
          len_nt <- (i - start_pos) * 3
          if (len_nt > best_orf$length) {
            start_genomic <- 1 + f + (start_pos - 1) * 3
            frame_val <- paste0("+", (start_genomic - 1) %% 3 + 1)
            
            best_orf <- list(
              length = len_nt,
              frame = frame_val,
              strand = "+",
              start = start_genomic,
              end = 1 + f + i * 3 - 1
            )
          }
          in_orf <- FALSE
        }
      }
    }
  }
  
  rev_seq <- reverse_complement(seq_upper)
  for (f in 0:2) {
    sub_seq <- rev_seq[(1 + f):L]
    n_codons <- floor(length(sub_seq) / 3)
    if (n_codons == 0) next
    
    idx_start <- seq(1, n_codons * 3, by = 3)
    codons <- paste0(sub_seq[idx_start], sub_seq[idx_start + 1], sub_seq[idx_start + 2])
    
    in_orf <- FALSE
    start_pos <- 0
    
    for (i in seq_along(codons)) {
      codon <- codons[i]
      if (!in_orf) {
        if (codon %in% starts) {
          in_orf <- TRUE
          start_pos <- i
        }
      } else {
        if (codon %in% stops) {
          len_nt <- (i - start_pos) * 3
          if (len_nt > best_orf$length) {
            rev_start_pos <- 1 + f + (start_pos - 1) * 3
            orig_5prime <- L - rev_start_pos + 1
            frame_num <- - ((orig_5prime %% 3) + 1)
            
            best_orf <- list(
              length = len_nt,
              frame = as.character(frame_num),
              strand = "-",
              rev_start = rev_start_pos,
              orig_5prime = orig_5prime
            )
          }
          in_orf <- FALSE
        }
      }
    }
  }
  
  return(best_orf)
}

sequence <- read.GenBank(access.nb = "NZ_JFIS01000023.1", as.character = TRUE)[[1]]

longueur_bp <- length(sequence)
real_orfs_lengths <- find_all_orfs(sequence, bacterial_starts, standard_stops)

set.seed(42)
null_sequence <- sample(sequence)
null_orfs_lengths <- find_all_orfs(null_sequence, bacterial_starts, standard_stops)

max_null_length <- max(null_orfs_lengths)
nb_candidate_genes <- sum(real_orfs_lengths > max_null_length)

threshold_p01 <- quantile(null_orfs_lengths, 0.99)
nb_candidates_p01 <- sum(real_orfs_lengths > threshold_p01)

longest_orf_res <- find_longest_candidate_frame(sequence)

print(longueur_bp)
print(max_null_length)
print(nb_candidate_genes)
print(nb_candidates_p01)
print(longest_orf_res)