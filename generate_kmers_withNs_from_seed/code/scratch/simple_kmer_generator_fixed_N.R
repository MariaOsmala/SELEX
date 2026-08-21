#!/usr/bin/env Rscript

# Simplified standalone script for generating k-mers with fixed N positions
# and Hamming distance 1 and 2 neighbors

# Function to expand IUPAC codes (N stays as N)
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
    'N' = c('N')  # N stays as N
  )
}

# Expand seed to all realizations (N positions remain as N)
expand_seed <- function(seed) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  iupac_map <- iupac_to_bases()
  
  sequences <- list("")
  
  for (char in seed_chars) {
    new_sequences <- list()
    
    if (char == 'N') {
      bases <- c('N')
    } else if (char %in% names(iupac_map)) {
      bases <- iupac_map[[char]]
    } else {
      stop(paste("Unknown character:", char))
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

# Generate Hamming-1 neighbors (only substitute at non-N positions)
generate_hamming1 <- function(sequence) {
  seq_chars <- strsplit(sequence, "")[[1]]
  neighbors <- c()
  bases <- c('A', 'C', 'G', 'T')
  
  for (i in seq_along(seq_chars)) {
    if (seq_chars[i] == 'N') next  # Skip N positions
    
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

# Generate Hamming-2 neighbors (only substitute at non-N positions)
generate_hamming2 <- function(sequence) {
  seq_chars <- strsplit(sequence, "")[[1]]
  neighbors <- c()
  bases <- c('A', 'C', 'G', 'T')
  
  non_n_pos <- which(seq_chars != 'N')
  
  if (length(non_n_pos) >= 2) {
    for (i in 1:(length(non_n_pos) - 1)) {
      for (j in (i + 1):length(non_n_pos)) {
        pos1 <- non_n_pos[i]
        pos2 <- non_n_pos[j]
        
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

# Main processing function
process_seed <- function(seed, output_prefix = "output") {
  cat("Processing seed:", seed, "\n")
  cat("=====================================\n")
  
  # Count positions
  seed_chars <- strsplit(toupper(seed), "")[[1]]
  n_positions <- which(seed_chars == 'N')
  non_n_positions <- which(seed_chars != 'N')
  
  cat("Seed length:", length(seed_chars), "\n")
  cat("N positions (fixed):", length(n_positions), "at positions:", paste(n_positions, collapse=", "), "\n")
  cat("Mutable positions:", length(non_n_positions), "at positions:", paste(non_n_positions, collapse=", "), "\n\n")
  
  # Step 1: Generate all realizations
  cat("Step 1: Generating realizations (N stays as N)...\n")
  realizations <- expand_seed(seed)
  cat("  Found", length(realizations), "realizations\n")
  
  # Show first few realizations as example
  if (length(realizations) <= 10) {
    cat("  All realizations:\n")
    for (r in realizations) {
      cat("    ", r, "\n")
    }
  } else {
    cat("  First 5 realizations:\n")
    for (i in 1:5) {
      cat("    ", realizations[i], "\n")
    }
    cat("  ... and", length(realizations) - 5, "more\n")
  }
  
  # Save realizations
  writeLines(realizations, paste0(output_prefix, "_realizations.txt"))
  cat("  Saved to:", paste0(output_prefix, "_realizations.txt"), "\n\n")
  
  # Step 2: Generate Hamming-1 neighbors
  cat("Step 2: Generating Hamming-1 neighbors...\n")
  all_hamming1 <- c()
  
  for (real in realizations) {
    h1 <- generate_hamming1(real)
    all_hamming1 <- c(all_hamming1, h1)
  }
  
  # Remove duplicates and exclude realizations
  all_hamming1 <- unique(all_hamming1)
  all_hamming1 <- setdiff(all_hamming1, realizations)
  
  cat("  Found", length(all_hamming1), "unique Hamming-1 neighbors\n")
  
  # Show examples
  if (length(all_hamming1) > 0) {
    cat("  First few Hamming-1 neighbors:\n")
    for (i in 1:min(5, length(all_hamming1))) {
      cat("    ", all_hamming1[i], "\n")
    }
  }
  
  # Save Hamming-1
  writeLines(all_hamming1, paste0(output_prefix, "_hamming1.txt"))
  cat("  Saved to:", paste0(output_prefix, "_hamming1.txt"), "\n\n")
  
  # Step 3: Generate Hamming-2 neighbors
  cat("Step 3: Generating Hamming-2 neighbors...\n")
  all_hamming2 <- c()
  
  for (real in realizations) {
    h2 <- generate_hamming2(real)
    all_hamming2 <- c(all_hamming2, h2)
  }
  
  # Remove duplicates and exclude realizations and Hamming-1
  all_hamming2 <- unique(all_hamming2)
  all_hamming2 <- setdiff(all_hamming2, c(realizations, all_hamming1))
  
  cat("  Found", length(all_hamming2), "unique Hamming-2 neighbors\n")
  
  # Show examples
  if (length(all_hamming2) > 0) {
    cat("  First few Hamming-2 neighbors:\n")
    for (i in 1:min(5, length(all_hamming2))) {
      cat("    ", all_hamming2[i], "\n")
    }
  }
  
  # Save Hamming-2
  writeLines(all_hamming2, paste0(output_prefix, "_hamming2.txt"))
  cat("  Saved to:", paste0(output_prefix, "_hamming2.txt"), "\n\n")
  
  # Summary
  cat("SUMMARY\n")
  cat("-------\n")
  cat("Seed:", seed, "\n")
  cat("Realizations:", length(realizations), "\n")
  cat("Hamming-1 neighbors:", length(all_hamming1), "(excluding realizations)\n")
  cat("Hamming-2 neighbors:", length(all_hamming2), "(excluding realizations and Hamming-1)\n")
  cat("Total unique sequences:", length(realizations) + length(all_hamming1) + length(all_hamming2), "\n")
  
  return(list(
    seed = seed,
    realizations = realizations,
    hamming1 = all_hamming1,
    hamming2 = all_hamming2
  ))
}

# Example usage
main <- function() {
  cat("K-mer Generation with Fixed N Positions\n")
  cat("========================================\n\n")
  
  # Example 1: Your provided seed
  cat("Example 1: Complex seed with multiple Ns\n")
  cat("-----------------------------------------\n")
  seed1 <- "NNCCGGNNNNNNCCGGNN"
  result1 <- process_seed(seed1, "example1")
  
  cat("\n\n")
  
  # Example 2: Simpler seed for demonstration
  cat("Example 2: Simple seed with IUPAC codes\n")
  cat("----------------------------------------\n")
  seed2 <- "ATNGCN"
  result2 <- process_seed(seed2, "example2")
  
  cat("\n\n")
  
  # Example 3: Seed with only IUPAC ambiguities
  cat("Example 3: Seed with IUPAC ambiguities\n")
  cat("---------------------------------------\n")
  seed3 <- "RRYYWS"
  result3 <- process_seed(seed3, "example3")
  
  return(list(
    example1 = result1,
    example2 = result2,
    example3 = result3
  ))
}

# Run the examples
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    # Process user-provided seed
    seed <- args[1]
    output_prefix <- ifelse(length(args) > 1, args[2], "output")
    result <- process_seed(seed, output_prefix)
  } else {
    # Run examples
    results <- main()
  }
} else {
  # If running interactively, show examples
  results <- main()
}
