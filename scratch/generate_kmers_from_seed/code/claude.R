#!/usr/bin/env Rscript

# Function to expand IUPAC codes to base possibilities
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
    'N' = c('A', 'C', 'G', 'T')
  )
}

# Function to recursively generate and stream realizations to file
stream_realizations <- function(seed_chars, position, current_seq, file_conn, iupac_map) {
  if (position > length(seed_chars)) {
    # Write the complete sequence to file (one per line)
    cat(paste(current_seq, collapse = ""), "\n", file = file_conn, sep = "")
    return(1)  # Return count of sequences written
  }
  
  char <- seed_chars[position]
  possible_bases <- if (char %in% names(iupac_map)) {
    iupac_map[[char]]
  } else {
    warning(paste("Unknown character:", char, "- treating as N"))
    c('A', 'C', 'G', 'T')
  }
  
  count <- 0
  for (base in possible_bases) {
    current_seq[position] <- base
    count <- count + stream_realizations(seed_chars, position + 1, current_seq, file_conn, iupac_map)
  }
  
  return(count)
}

# Function to stream Hamming-1 neighbors for a given sequence
stream_hamming1_for_sequence <- function(sequence, seed_chars, non_n_positions, file_conn, written_seqs = NULL) {
  bases <- c('A', 'C', 'G', 'T')
  seq_chars <- strsplit(sequence, "")[[1]]
  count <- 0
  
  # For each non-N position, try all possible substitutions
  for (pos in non_n_positions) {
    current_base <- seq_chars[pos]
    
    # Try each alternative base
    for (new_base in bases) {
      if (new_base != current_base) {
        # Create new sequence with substitution
        new_seq <- seq_chars
        new_seq[pos] <- new_base
        new_seq_str <- paste(new_seq, collapse = "")
        
        # Check if we've already written this sequence (for deduplication)
        if (is.null(written_seqs) || !new_seq_str %in% written_seqs) {
          cat(new_seq_str, "\n", file = file_conn, sep = "")
          count <- count + 1
          if (!is.null(written_seqs)) {
            written_seqs <- c(written_seqs, new_seq_str)
          }
        }
      }
    }
  }
  
  return(list(count = count, written = written_seqs))
}

# Function to process realizations in chunks and generate Hamming-1 neighbors
process_realizations_chunked <- function(seed_chars, non_n_positions, realizations_file, hamming_file, 
                                         chunk_size = 10000, verbose = TRUE) {
  iupac_map <- iupac_to_bases()
  
  # Read realizations in chunks and process
  conn_in <- file(realizations_file, "r")
  conn_out <- file(hamming_file, "w")
  
  written_hamming <- character()  # Track written sequences for deduplication
  total_hamming <- 0
  total_realizations <- 0
  
  repeat {
    # Read chunk of realizations
    lines <- readLines(conn_in, n = chunk_size)
    if (length(lines) == 0) break
    
    total_realizations <- total_realizations + length(lines)
    
    # Process each realization in the chunk
    for (realization in lines) {
      result <- stream_hamming1_for_sequence(realization, seed_chars, non_n_positions, 
                                             conn_out, written_hamming)
      written_hamming <- result$written
      total_hamming <- total_hamming + result$count
    }
    
    if (verbose && total_realizations %% 10000 == 0) {
      cat("  Processed", total_realizations, "realizations,", 
          "generated", total_hamming, "Hamming-1 neighbors\n")
    }
  }
  
  close(conn_in)
  close(conn_out)
  
  return(list(realizations = total_realizations, hamming1 = total_hamming))
}

# Main function to process a seed sequence with streaming
process_seed_streaming <- function(seed, output_prefix = NULL, verbose = TRUE, 
                                   deduplicate_hamming = TRUE) {
  if (verbose) {
    cat("Processing seed:", seed, "\n")
    cat("========================================\n")
  }
  
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  iupac_map <- iupac_to_bases()
  
  # Find non-N positions for Hamming distance calculation
  non_n_positions <- which(seed_chars != 'N')
  
  if (verbose) {
    cat("Seed length:", length(seed_chars), "\n")
    cat("Non-N positions:", length(non_n_positions), "\n")
    cat("Positions that can be mutated:", paste(non_n_positions, collapse = ", "), "\n\n")
  }
  
  # Generate output filenames
  if (is.null(output_prefix)) {
    output_prefix <- gsub("[^A-Za-z0-9]", "_", seed)
  }
  
  realizations_file <- paste0(output_prefix, "_realizations.txt")
  hamming_file <- paste0(output_prefix, "_hamming1.txt")
  
  # Step 1: Stream all realizations to file
  if (verbose) {
    cat("Streaming realizations to:", realizations_file, "\n")
  }
  
  conn_real <- file(realizations_file, "w")
  realization_count <- stream_realizations(seed_chars, 1, character(length(seed_chars)), 
                                           conn_real, iupac_map)
  close(conn_real)
  
  if (verbose) {
    cat("  Total realizations written:", realization_count, "\n\n")
  }
  
  # Step 2: Stream Hamming-1 neighbors
  if (verbose) {
    cat("Generating Hamming-1 neighbors to:", hamming_file, "\n")
  }
  
  if (deduplicate_hamming) {
    # Use hash table for deduplication (more memory but no duplicates)
    hamming_count <- generate_hamming1_deduplicated(realizations_file, hamming_file, 
                                                    seed_chars, non_n_positions, verbose)
  } else {
    # Simple streaming (less memory but may have duplicates)
    counts <- process_realizations_chunked(seed_chars, non_n_positions, 
                                           realizations_file, hamming_file, verbose = verbose)
    hamming_count <- counts$hamming1
  }
  
  if (verbose) {
    cat("  Total Hamming-1 neighbors written:", hamming_count, "\n\n")
  }
  
  return(list(
    seed = seed,
    realizations_file = realizations_file,
    hamming1_file = hamming_file,
    realization_count = realization_count,
    hamming1_count = hamming_count
  ))
}

# Function to generate deduplicated Hamming-1 neighbors using a hash table
generate_hamming1_deduplicated <- function(realizations_file, hamming_file, 
                                           seed_chars, non_n_positions, verbose = TRUE) {
  bases <- c('A', 'C', 'G', 'T')
  seen_sequences <- new.env(hash = TRUE)  # Use environment as hash table
  
  conn_in <- file(realizations_file, "r")
  conn_out <- file(hamming_file, "w")
  
  total_written <- 0
  total_processed <- 0
  chunk_size <- 10000
  
  repeat {
    lines <- readLines(conn_in, n = chunk_size)
    if (length(lines) == 0) break
    
    for (realization in lines) {
      seq_chars <- strsplit(realization, "")[[1]]
      
      # Generate all Hamming-1 neighbors
      for (pos in non_n_positions) {
        current_base <- seq_chars[pos]
        
        for (new_base in bases) {
          if (new_base != current_base) {
            new_seq <- seq_chars
            new_seq[pos] <- new_base
            new_seq_str <- paste(new_seq, collapse = "")
            
            # Check if we've seen this sequence before
            if (is.null(seen_sequences[[new_seq_str]])) {
              seen_sequences[[new_seq_str]] <- TRUE
              cat(new_seq_str, "\n", file = conn_out, sep = "")
              total_written <- total_written + 1
            }
          }
        }
      }
      
      total_processed <- total_processed + 1
      if (verbose && total_processed %% 10000 == 0) {
        cat("  Processed", total_processed, "realizations,",
            "unique Hamming-1 neighbors:", total_written, "\n")
      }
    }
  }
  
  close(conn_in)
  close(conn_out)
  
  return(total_written)
}

# Function to process multiple seeds with streaming
process_multiple_seeds_streaming <- function(seeds, output_dir = ".", verbose = TRUE,
                                             deduplicate_hamming = TRUE) {
  results <- list()
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (i in seq_along(seeds)) {
    if (verbose) {
      cat("\n=== Processing seed", i, "of", length(seeds), "===\n")
    }
    
    # Generate output prefix with path
    output_prefix <- file.path(output_dir, paste0("seed_", i))
    
    results[[i]] <- process_seed_streaming(seeds[i], output_prefix, verbose, deduplicate_hamming)
  }
  
  names(results) <- seeds
  return(results)
}






# Example usage and main function
library("readr")
library("dplyr")
metadata=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/metadata_failed_and_successfull_seeds_degeneracy.tsv")

motif=metadata %>% filter(kmer_generation_successfull==TRUE) %>% filter(length==(min(length)+1)) %>% pull(ID)
motif=motif[1]
s=metadata %>% filter(ID==motif) %>% pull(seed)

# Example with the provided seed

seed1 <- s
seed=seed1
result1 <- process_seed_streaming(seed1, output_prefix=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif), 
                                    verbose = TRUE, deduplicate_hamming = TRUE)

s=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(seed)
motif=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(ID)
#print(seq_stats(s))
result1 <- process_seed_streaming(s, output_prefix=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif), 
                                  verbose = TRUE, deduplicate_hamming = TRUE)
#List of 5
#$ k                : int 18
#$ mutable_positions: int [1:8] 3 4 5 6 13 14 15 16
#$ n_realizations   : num 1048576
#$ hamming1_per_real: int 24
#$ total_hamming1   : num 25,165,824
  
  # Example with multiple seeds
  cat("\n### EXAMPLE 2: Multiple seeds with streaming ###\n")
  seeds <- c(
    "NNCCGGNNNNNNCCGGNN",
    "ATNGCN",
    "RRYYWS"
  )
  
  all_results <- process_multiple_seeds_streaming(
    seeds, 
    output_dir = "kmer_output",
    verbose = TRUE,
    deduplicate_hamming = TRUE
  )
  
  # Summary
  cat("\n### SUMMARY ###\n")
  cat("----------------------------------------\n")
  for (seed in names(all_results)) {
    cat("Seed:", seed, "\n")
    cat("  Realizations file:", all_results[[seed]]$realizations_file, "\n")
    cat("  Realizations count:", all_results[[seed]]$realization_count, "\n")
    cat("  Hamming-1 file:", all_results[[seed]]$hamming1_file, "\n")
    cat("  Hamming-1 count:", all_results[[seed]]$hamming1_count, "\n")
    cat("----------------------------------------\n")
  }
  
  return(all_results)


# Function to process seeds from an input file
process_seeds_from_file <- function(input_file, output_dir = "output", 
                                    verbose = TRUE, deduplicate_hamming = TRUE) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist:", input_file)
  }
  
  # Read seeds from file
  seeds <- readLines(input_file)
  seeds <- seeds[seeds != ""]  # Remove empty lines
  
  cat("Found", length(seeds), "seeds in", input_file, "\n\n")
  
  # Process all seeds
  results <- process_multiple_seeds_streaming(seeds, output_dir, verbose, deduplicate_hamming)
  
  # Write summary file
  summary_file <- file.path(output_dir, "processing_summary.txt")
  sink(summary_file)
  cat("Processing Summary\n")
  cat("==================\n")
  cat("Date:", date(), "\n")
  cat("Input file:", input_file, "\n")
  cat("Total seeds processed:", length(seeds), "\n\n")
  
  for (seed in names(results)) {
    cat("Seed:", seed, "\n")
    cat("  Realizations:", results[[seed]]$realization_count, "\n")
    cat("  Hamming-1 neighbors:", results[[seed]]$hamming1_count, "\n")
    cat("\n")
  }
  sink()
  
  cat("\nSummary written to:", summary_file, "\n")
  
  return(results)
}

# Run if script is executed directly
if (!interactive()) {
  # Check command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    # Process seeds from input file
    input_file <- args[1]
    output_dir <- ifelse(length(args) > 1, args[2], "output")
    results <- process_seeds_from_file(input_file, output_dir)
  } else {
    # Run example
    results <- main()
  }
}