#!/usr/bin/env Rscript

# Function to expand IUPAC codes to base possibilities (N stays as N)
iupac_to_bases <- function() {
  list(
    'A' = c('A'),
    'C' = c('C'),
    'G' = c('G'),
    'T' = c('T'),
    'R' = c('A', 'G'),
    'Y' = c('C', 'T'),
    'S' = c('G', 'C'),
    'W' = c('A', 'T'),
    'K' = c('G', 'T'),
    'M' = c('A', 'C'),
    'B' = c('C', 'G', 'T'),
    'D' = c('A', 'G', 'T'),
    'H' = c('A', 'C', 'T'),
    'V' = c('A', 'C', 'G'),
    'N' = c('N')  # CHANGED: N stays as N instead of expanding to all 4 bases
  )
}

# Function to compute reverse complement (N stays as N)
reverse_complement <- function(sequence) {
  complement_map <- c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C', 'N' = 'N')
  seq_chars <- strsplit(sequence, "")[[1]]
  rev_comp_chars <- complement_map[rev(seq_chars)]
  return(paste(rev_comp_chars, collapse = ""))
}

# Function to get canonical k-mer (lexicographically smallest)
get_canonical <- function(sequence) {
  rev_comp <- reverse_complement(sequence)
  if (sequence <= rev_comp) {
    return(sequence)
  } else {
    return(rev_comp)
  }
}



# Function to expand seed to all realizations (with N fixed)
expand_seed <- function(seed) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  iupac_map <- iupac_to_bases()
  
  # Start with a list containing an empty sequence
  sequences <- list("")
  
  for (char in seed_chars) {
    new_sequences <- list()
    
    if (char == 'N') {
      # N stays as N
      bases <- c('N')
    } else if (char %in% names(iupac_map)) {
      bases <- iupac_map[[char]]
    } else {
      warning(paste("Unknown character:", char, "- skipping"))
      bases <- c()
    }
    
    for (seq in sequences) {
      for (base in bases) {
        new_sequences <- append(new_sequences, paste0(seq, base))
      }
    }
    
    sequences <- new_sequences
  }
  
  return(unlist(sequences))
}

# Function to generate Hamming distance 1 neighbors (substitutions only at non-N positions)
generate_hamming1_neighbors <- function(sequence) {
  seq_chars <- strsplit(sequence, "")[[1]]
  neighbors <- c()
  bases <- c('A', 'C', 'G', 'T')
  
  for (i in seq_along(seq_chars)) {
    # Skip N positions - no substitutions allowed here
    if (seq_chars[i] == 'N') {
      next
    }
    
    # Generate all possible single substitutions at this position
    for (base in bases) {
      if (base != seq_chars[i]) {
        new_seq <- seq_chars
        new_seq[i] <- base
        neighbors <- c(neighbors, paste(new_seq, collapse = ""))
      }
    }
  }
  
  return(unique(neighbors))
}

# Function to generate Hamming distance 2 neighbors (substitutions only at non-N positions)
generate_hamming2_neighbors <- function(sequence) {
  seq_chars <- strsplit(sequence, "")[[1]]
  neighbors <- c()
  bases <- c('A', 'C', 'G', 'T')
  
  # Find all non-N positions
  non_n_positions <- which(seq_chars != 'N')
  
  # Generate all pairs of positions for double substitutions
  if (length(non_n_positions) >= 2) {
    for (i in 1:(length(non_n_positions) - 1)) {
      for (j in (i + 1):length(non_n_positions)) {
        pos1 <- non_n_positions[i]
        pos2 <- non_n_positions[j]
        
        # Try all combinations of bases at these two positions
        for (base1 in bases) {
          if (base1 != seq_chars[pos1]) {
            for (base2 in bases) {
              if (base2 != seq_chars[pos2]) {
                new_seq <- seq_chars
                new_seq[pos1] <- base1
                new_seq[pos2] <- base2
                neighbors <- c(neighbors, paste(new_seq, collapse = ""))
              }
            }
          }
        }
      }
    }
  }
  
  return(unique(neighbors))
}

# Streaming function to process a single seed
process_seed_streaming <- function(seed, output_prefix = "output", 
                                  verbose = TRUE, 
                                  canonical_only = FALSE,
                                  generate_hamming2 = TRUE) {
  
  seed_upper <- toupper(seed)
  
  if (verbose) {
    cat("Processing seed:", seed_upper, "\n")
    cat("Seed length:", nchar(seed_upper), "\n")
    
    # Count N positions
    n_count <- sum(strsplit(seed_upper, "")[[1]] == 'N')
    cat("Number of N positions (fixed):", n_count, "\n")
    cat("Number of mutable positions:", nchar(seed_upper) - n_count, "\n")
  }
  
  # Generate output file names
  if (canonical_only) {
    realizations_file <- paste0(output_prefix, "_canonical_realizations.txt")
    hamming1_file <- paste0(output_prefix, "_canonical_hamming1.txt")
    hamming2_file <- paste0(output_prefix, "_canonical_hamming2.txt")
  } else {
    realizations_file <- paste0(output_prefix, "_realizations.txt")
    hamming1_file <- paste0(output_prefix, "_hamming1.txt")
    hamming2_file <- paste0(output_prefix, "_hamming2.txt")
  }
  
  # Expand seed to get all realizations (with N fixed)
  if (verbose) cat("\nExpanding seed to realizations (N positions stay as N)...\n")
  realizations <- expand_seed(seed_upper)
  
  if (canonical_only) {
    realizations <- unique(sapply(realizations, get_canonical))
  }
  
  if (verbose) cat("  Generated", length(realizations), "realizations\n")
  
  # Write realizations
  writeLines(realizations, realizations_file)
  if (verbose) cat("  Written to:", realizations_file, "\n")
  
  # Track all unique sequences for deduplication
  all_hamming1 <- c()
  all_hamming2 <- c()
  realizations_set <- as.list(setNames(rep(TRUE, length(realizations)), realizations))
  
  if (verbose) cat("\nGenerating Hamming distance 1 neighbors...\n")
  
  # Generate Hamming-1 neighbors for each realization
  for (realization in realizations) {
    h1_neighbors <- generate_hamming1_neighbors(realization)
    
    if (canonical_only) {
      h1_neighbors <- sapply(h1_neighbors, get_canonical)
    }
    
    # Exclude sequences that are realizations
    h1_neighbors <- h1_neighbors[!(h1_neighbors %in% names(realizations_set))]
    
    all_hamming1 <- c(all_hamming1, h1_neighbors)
  }
  
  # Deduplicate Hamming-1 neighbors
  all_hamming1 <- unique(all_hamming1)
  
  # Write Hamming-1 neighbors
  if (length(all_hamming1) > 0) {
    writeLines(all_hamming1, hamming1_file)
    if (verbose) cat("  Generated", length(all_hamming1), "unique Hamming-1 neighbors\n")
    if (verbose) cat("  Written to:", hamming1_file, "\n")
  } else {
    if (verbose) cat("  No Hamming-1 neighbors found\n")
  }
  
  # Generate Hamming-2 neighbors if requested
  if (generate_hamming2) {
    if (verbose) cat("\nGenerating Hamming distance 2 neighbors...\n")
    
    # Create a set of sequences to exclude (realizations + Hamming-1)
    exclude_set <- c(names(realizations_set), all_hamming1)
    exclude_set <- as.list(setNames(rep(TRUE, length(exclude_set)), exclude_set))
    
    for (realization in realizations) {
      h2_neighbors <- generate_hamming2_neighbors(realization)
      
      if (canonical_only) {
        h2_neighbors <- sapply(h2_neighbors, get_canonical)
      }
      
      # Exclude sequences that are realizations or Hamming-1 neighbors
      h2_neighbors <- h2_neighbors[!(h2_neighbors %in% names(exclude_set))]
      
      all_hamming2 <- c(all_hamming2, h2_neighbors)
    }
    
    # Deduplicate Hamming-2 neighbors
    all_hamming2 <- unique(all_hamming2)
    
    # Write Hamming-2 neighbors
    if (length(all_hamming2) > 0) {
      writeLines(all_hamming2, hamming2_file)
      if (verbose) cat("  Generated", length(all_hamming2), "unique Hamming-2 neighbors\n")
      if (verbose) cat("  Written to:", hamming2_file, "\n")
    } else {
      if (verbose) cat("  No Hamming-2 neighbors found\n")
    }
  }
  
  # Return summary
  result <- list(
    seed = seed_upper,
    canonical_only = canonical_only,
    realizations_file = realizations_file,
    realization_count = length(realizations),
    hamming1_file = hamming1_file,
    hamming1_count = length(all_hamming1)
  )
  
  if (generate_hamming2) {
    result$hamming2_file <- hamming2_file
    result$hamming2_count <- length(all_hamming2)
  }
  
  return(result)
}







# Modified version for HPC cluster usage (if using with your metadata)

library("readr")
library("dplyr")
  
args <- commandArgs(trailingOnly = TRUE)
  
arrays <- as.numeric(args[1])
addition <- as.numeric(args[2])
arrays <- arrays + addition
  
# Load metadata
metadata <- read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim = "\t") #3933
metadata <- metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata <- metadata %>% filter(!is.na(seed)) #3635
  
# Sort based on seed length
metadata <- metadata[order(metadata$length),]
  
length <- 10
start_ind <- (arrays - 1) * length + 1
end_ind <- arrays * length
  
if (end_ind > nrow(metadata)) {
  end_ind <- nrow(metadata)
}
  
path_local_scratch <- Sys.getenv("LOCAL_SCRATCH")
lustre_folder <- "/scratch/project_2013895/SELEX/kmers_from_fixed_N_seeds/"
  
for (i in start_ind:end_ind) {
  
  print(i)
  motif <- metadata$ID[i]
  seed <- metadata$seed[i]
  
  #which(metadata$seed=="NNCCGGNNNNNNCCGGNN")#
  #i=3049
  
  result <- process_seed_streaming(
    seed, 
    output_prefix = paste0(path_local_scratch, "/", motif),
    verbose = TRUE,
    canonical_only = TRUE,
    generate_hamming2 = TRUE
  )
  
  print(str(result))
  
  # Copy files to lustre
  system(paste0("cp ", path_local_scratch, "/", motif, "_canonical_realizations.txt ", lustre_folder))
  system(paste0("cp ", path_local_scratch, "/", motif, "_canonical_hamming1.txt ", lustre_folder))
  system(paste0("cp ", path_local_scratch, "/", motif, "_canonical_hamming2.txt ", lustre_folder))
}



