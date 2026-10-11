library(seqinr) # ou library(Biostrings) pour la lecture FASTA

bacterial_starts <- c("TTG", "CTG", "ATA", "ATT", "ATC", "ATG", "GTG")
standard_stops   <- c("TGA", "TAA", "TAG")
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
  
  all_orfs <- c(fwd_orfs, rev_orfs)
  return(all_orfs)
}

read_fasta_native <- function(file_path) {
  lines <- readLines(file_path)
  seq_lines <- lines[!grepl("^>", lines)]
  full_seq <- paste0(seq_lines, collapse = "")
  return(strsplit(full_seq, "")[[1]])
}

analyze_bacterial_genome <- function(fasta_file, output_rds = "orf_counts.rds") {
  sequence <- read_fasta_native(fasta_file)
  orfs_lengths <- find_all_orfs(sequence, bacterial_starts, standard_stops)
  
  k_values <- c(0, 10, 50, 100, 300, 500)  
  counts <- sapply(k_values, function(k) sum(orfs_lengths >= k))
  names(counts) <- as.character(k_values)  
  saveRDS(counts, file = output_rds)
  
  cat("Analyse terminée. Vecteur sauvegardé dans :", output_rds, "\n")
  return(counts)
}