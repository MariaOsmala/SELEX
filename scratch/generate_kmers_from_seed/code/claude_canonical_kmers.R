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

# Function to compute reverse complement
reverse_complement <- function(sequence) {
  complement_map <- c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C')
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

# ============================================================================
# THEORETICAL COUNTING FUNCTIONS
# ============================================================================

# Function to count total possible realizations (non-canonical)
count_realizations_theoretical <- function(seed) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  iupac_map <- iupac_to_bases()
  
  count <- 1
  for (char in seed_chars) {
    if (char %in% names(iupac_map)) {
      count <- count * length(iupac_map[[char]])
    } else {
      warning(paste("Unknown character:", char, "- treating as N"))
      count <- count * 4
    }
  }
  
  return(count)
}

# Function to estimate canonical realizations (accounting for palindromes)
count_canonical_realizations_theoretical <- function(seed) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  iupac_map <- iupac_to_bases()
  
  # For exact count, we need to consider palindromes
  # This is a simplified estimate that assumes few palindromes
  total_realizations <- count_realizations_theoretical(seed)
  
  # Check if the seed pattern could produce palindromes
  seed_length <- length(seed_chars)
  
  if (seed_length %% 2 == 1) {
    # Odd length - palindromes are possible but rare
    # Rough estimate: divide by 2 and add back a small fraction for palindromes
    canonical_estimate <- ceiling(total_realizations / 2)
  } else {
    # Even length - palindromes are possible
    # For a more accurate estimate, we'd need to check the seed pattern
    # This is still an approximation
    canonical_estimate <- ceiling(total_realizations / 2)
  }
  
  return(list(
    total_realizations = total_realizations,
    canonical_estimate = canonical_estimate,
    note = "Canonical count is an estimate; exact count requires checking for palindromic sequences"
  ))
}

# Function to count theoretical Hamming-1 neighbors
count_hamming1_theoretical <- function(seed, exclude_realizations = TRUE) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  
  # Find non-N positions
  non_n_positions <- which(seed_chars != 'N')
  num_mutable_positions <- length(non_n_positions)
  
  # Get number of realizations
  num_realizations <- count_realizations_theoretical(seed)
  
  # Each realization can have up to 3 * num_mutable_positions Hamming-1 neighbors
  # (3 alternative bases for each mutable position)
  max_hamming1_per_realization <- 3 * num_mutable_positions
  
  # Total theoretical Hamming-1 neighbors (with possible overlaps)
  total_hamming1 <- num_realizations * max_hamming1_per_realization
  
  # Calculate overlap probability
  overlap_info <- calculate_hamming1_overlap_probability(seed, non_n_positions)
  
  return(list(
    num_realizations = num_realizations,
    mutable_positions = num_mutable_positions,
    max_hamming1_per_realization = max_hamming1_per_realization,
    total_hamming1_with_overlaps = total_hamming1,
    overlap_probability = overlap_info$overlap_prob,
    estimated_unique_hamming1 = overlap_info$estimated_unique,
    note = if(exclude_realizations) 
      "If excluding realizations from Hamming-1, actual count will be lower" 
    else 
      "Includes potential overlaps between Hamming-1 neighbors"
  ))
}

# Function to calculate overlap probability for Hamming-1 neighbors
calculate_hamming1_overlap_probability <- function(seed, non_n_positions) {
  seed_chars <- strsplit(toupper(seed), "")[[1]]
  iupac_map <- iupac_to_bases()
  
  # Count variable positions (N or ambiguous IUPAC)
  variable_positions <- sapply(seed_chars, function(char) {
    if (char %in% names(iupac_map)) {
      length(iupac_map[[char]])
    } else {
      4
    }
  })
  
  num_realizations <- prod(variable_positions)
  num_mutable <- length(non_n_positions)
  
  # Probability that two random Hamming-1 neighbors are the same
  # This happens when two different mutations lead to the same sequence
  # Simplified calculation - this is an approximation
  
  # If we have many N positions, overlap is more likely
  n_positions <- which(seed_chars == 'N')
  num_n <- length(n_positions)
  
  # Overlaps can occur when:
  # 1. Different realizations have Hamming-1 neighbors that converge
  # 2. A Hamming-1 neighbor of one realization equals another realization
  
  # Rough estimate of unique Hamming-1 neighbors
  if (num_n == 0) {
    # No N positions - less overlap
    overlap_prob <- 0.05  # Rough estimate
    estimated_unique <- round(num_realizations * 3 * num_mutable * (1 - overlap_prob))
  } else {
    # With N positions, more overlap possible
    overlap_factor <- min(0.5, num_n * 0.1)  # Increases with more Ns
    overlap_prob <- overlap_factor
    estimated_unique <- round(num_realizations * 3 * num_mutable * (1 - overlap_prob))
  }
  
  return(list(
    overlap_prob = overlap_prob,
    estimated_unique = estimated_unique
  ))
}

# Function to get comprehensive theoretical counts
theoretical_counts <- function(seed, canonical = TRUE, verbose = TRUE) {
  seed_upper <- toupper(seed)
  seed_chars <- strsplit(seed_upper, "")[[1]]
  
  if (verbose) {
    cat("Theoretical counts for seed:", seed, "\n")
    cat("========================================\n")
    cat("Seed length:", length(seed_chars), "\n")
  }
  
  # Count each type of position
  # Count how many times each IUPAC character appears
  position_types <- table(seed_chars)
  
  if (verbose) {
    cat("\nPosition composition:\n")
    for (char in names(position_types)) {
      cat("  ", char, ":", position_types[char], "positions\n")
    }
  }
  
  # Get realization counts
  if (canonical) {
    real_counts <- count_canonical_realizations_theoretical(seed)
    realization_count <- real_counts$canonical_estimate
    
    if (verbose) {
      cat("\nRealizations:\n")
      cat("  Total (non-canonical):", real_counts$total_realizations, "\n")
      cat("  Canonical (estimate):", real_counts$canonical_estimate, "\n")
    }
  } else {
    realization_count <- count_realizations_theoretical(seed)
    
    if (verbose) {
      cat("\nRealizations:\n")
      cat("  Total:", realization_count, "\n")
    }
  }
  
  # Get Hamming-1 counts
  hamming_info <- count_hamming1_theoretical(seed)
  
  if (verbose) {
    cat("\nHamming-1 neighbors:\n")
    cat("  Mutable positions (non-N):", hamming_info$mutable_positions, "\n")
    cat("  Max per realization:", hamming_info$max_hamming1_per_realization, "\n")
    cat("  Total with overlaps:", hamming_info$total_hamming1_with_overlaps, "\n")
    cat("  Estimated unique:", hamming_info$estimated_unique_hamming1, "\n")
    
    if (canonical) {
      canonical_hamming_estimate <- ceiling(hamming_info$estimated_unique_hamming1 / 2)
      cat("  Canonical estimate:", canonical_hamming_estimate, "\n")
    }
  }
  
  return(list(
    seed = seed,
    seed_length = length(seed_chars),
    position_types = position_types,
    canonical_mode = canonical,
    realization_count = realization_count,
    hamming1_info = hamming_info
  ))
}

# Function to compare theoretical vs actual counts
compare_theoretical_actual <- function(seed, verbose = TRUE) {
  cat("Comparing theoretical vs actual counts\n")
  cat("======================================\n\n")
  
  # Get theoretical counts
  cat("THEORETICAL COUNTS:\n")
  theory <- theoretical_counts(seed, canonical = TRUE, verbose = FALSE)
  cat("  Canonical realizations (estimate):", theory$realization_count, "\n")
  cat("  Unique Hamming-1 (estimate):", theory$hamming1_info$estimated_unique_hamming1, "\n")
  cat("  Canonical Hamming-1 (estimate):", ceiling(theory$hamming1_info$estimated_unique_hamming1 / 2), "\n\n")
  
  # Generate actual counts
  cat("ACTUAL COUNTS:\n")
  cat("Generating actual sequences...\n")
  
  # Create temporary files
  temp_prefix <- tempfile(pattern = "compare_")
  
  result <- process_seed_streaming(
    seed = seed,
    output_prefix = temp_prefix,
    verbose = FALSE,
    canonical_only = TRUE,
    exclude_realizations_from_hamming = TRUE
  )
  
  cat("  Canonical realizations (actual):", result$realization_count, "\n")
  cat("  Canonical Hamming-1 (actual):", result$hamming1_count, "\n\n")
  
  # Calculate accuracy
  cat("ACCURACY:\n")
  real_accuracy <- round(100 * min(theory$realization_count, result$realization_count) / 
                           max(theory$realization_count, result$realization_count), 2)
  cat("  Realization estimate accuracy:", real_accuracy, "%\n")
  
  hamming_estimate <- ceiling(theory$hamming1_info$estimated_unique_hamming1 / 2)
  hamming_accuracy <- round(100 * min(hamming_estimate, result$hamming1_count) / 
                              max(hamming_estimate, result$hamming1_count), 2)
  cat("  Hamming-1 estimate accuracy:", hamming_accuracy, "%\n")
  
  # Clean up temp files
  unlink(result$realizations_file)
  unlink(result$hamming1_file)
  
  return(invisible(list(
    theoretical = theory,
    actual = result
  )))
}

# Function to recursively generate and stream realizations to file with tracking
stream_realizations_with_tracking <- function(seed_chars, position, current_seq, file_conn, iupac_map, 
                                              canonical_only = FALSE, seen_canonical = NULL, track_set = NULL) {
  
  #seed_chars
  #position=1
  #current_seq=character(length(seed_chars)), 
  #file_conn=conn_real
  #canonical_only = TRUE 
  #track_set = realizations_set
  
  
  if (position > length(seed_chars)) {
    # Get the complete sequence
    sequence <- paste(current_seq, collapse = "")
    
    # If canonical_only, check and write only canonical form
    if (canonical_only) {
      canonical_seq <- get_canonical(sequence)
      # Check if we've already written this canonical sequence
      if (is.null(seen_canonical[[canonical_seq]])) {
        seen_canonical[[canonical_seq]] <- TRUE
        cat(canonical_seq, "\n", file = file_conn, sep = "")
        # Also track in the realizations set if needed
        if (!is.null(track_set)) {
          track_set[[canonical_seq]] <- TRUE
        }
        return(list(count = 1, seen = seen_canonical))
      } else {
        return(list(count = 0, seen = seen_canonical))
      }
    } else {
      # Write the sequence as-is
      cat(sequence, "\n", file = file_conn, sep = "")
      # Track in the realizations set if needed
      if (!is.null(track_set)) {
        track_set[[sequence]] <- TRUE
      }
      return(list(count = 1, seen = seen_canonical))
    }
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
    result <- stream_realizations_with_tracking(seed_chars, position + 1, current_seq, file_conn, 
                                                iupac_map, canonical_only, seen_canonical, track_set)
    count <- count + result$count
    seen_canonical <- result$seen
  }
  
  return(list(count = count, seen = seen_canonical))
}

# Function to stream Hamming-1 neighbors for a given sequence
stream_hamming1_for_sequence <- function(sequence, seed_chars, non_n_positions, file_conn, 
                                         written_seqs = NULL, canonical_only = FALSE,
                                         realizations_set = NULL) {
  bases <- c('A', 'C', 'G', 'T')
  seq_chars <- strsplit(sequence, "")[[1]]
  count <- 0
  excluded <- 0
  
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
        
        # Get canonical form if requested
        if (canonical_only) {
          new_seq_str <- get_canonical(new_seq_str)
        }
        
        # Check if this sequence is a realization (and should be excluded)
        if (!is.null(realizations_set) && !is.null(realizations_set[[new_seq_str]])) {
          excluded <- excluded + 1
          next  # Skip this sequence
        }
        
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
  
  return(list(count = count, written = written_seqs, excluded = excluded))
}

# Function to process realizations in chunks and generate Hamming-1 neighbors
process_realizations_chunked <- function(seed_chars, non_n_positions, realizations_file, hamming_file, 
                                         chunk_size = 10000, verbose = TRUE, canonical_only = FALSE,
                                         realizations_set = NULL) {
  iupac_map <- iupac_to_bases()
  
  # Read realizations in chunks and process
  conn_in <- file(realizations_file, "r")
  conn_out <- file(hamming_file, "w")
  
  written_hamming <- character()  # Track written sequences for deduplication
  total_hamming <- 0
  total_realizations <- 0
  total_excluded <- 0
  
  repeat {
    # Read chunk of realizations
    lines <- readLines(conn_in, n = chunk_size)
    if (length(lines) == 0) break
    
    total_realizations <- total_realizations + length(lines)
    
    # Process each realization in the chunk
    for (realization in lines) {
      result <- stream_hamming1_for_sequence(realization, seed_chars, non_n_positions, 
                                             conn_out, written_hamming, canonical_only,
                                             realizations_set)
      written_hamming <- result$written
      total_hamming <- total_hamming + result$count
      if (!is.null(result$excluded)) {
        total_excluded <- total_excluded + result$excluded
      }
    }
    
    if (verbose && total_realizations %% 10000 == 0) {
      cat("  Processed", total_realizations, "realizations,", 
          "generated", total_hamming, "Hamming-1 neighbors")
      if (!is.null(realizations_set)) {
        cat(" (excluded", total_excluded, ")")
      }
      cat("\n")
    }
  }
  
  close(conn_in)
  close(conn_out)
  
  return(list(realizations = total_realizations, hamming1 = total_hamming))
}

# Main function to process a seed sequence with streaming
process_seed_streaming <- function(seed, output_prefix = NULL, verbose = TRUE, 
                                   deduplicate_hamming = TRUE, canonical_only = FALSE,
                                   exclude_realizations_from_hamming = TRUE) {
  
  #seed=seed1
  #output_prefix=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif)
  #verbose = TRUE
  #deduplicate_hamming = TRUE
  #canonical_only=TRUE
  #exclude_realizations_from_hamming = TRUE
  
  
  
  if (verbose) {
    cat("Processing seed:", seed, "\n")
    if (canonical_only) {
      cat("Mode: Canonical k-mers only (lexicographically smallest)\n")
    }
    if (exclude_realizations_from_hamming) {
      cat("Excluding realizations from Hamming-1 output\n")
    }
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
  
  if (canonical_only) {
    realizations_file <- paste0(output_prefix, "_canonical_realizations.txt")
    hamming_file <- paste0(output_prefix, "_canonical_hamming1.txt")
  } else {
    realizations_file <- paste0(output_prefix, "_realizations.txt")
    hamming_file <- paste0(output_prefix, "_hamming1.txt")
  }
  
  # Step 1: Stream all realizations to file
  if (verbose) {
    cat("Streaming realizations to:", realizations_file, "\n")
  }
  
  conn_real <- file(realizations_file, "w")
  
  # Keep track of realizations if we need to exclude them from Hamming-1
  realizations_set <- NULL
  if (exclude_realizations_from_hamming) {
    realizations_set <- new.env(hash = TRUE)
  }
  
  if (canonical_only) {
    # Use hash table to track canonical sequences
    seen_canonical <- new.env(hash = TRUE)
    result <- stream_realizations_with_tracking(seed_chars, 1, character(length(seed_chars)), 
                                                conn_real, iupac_map, canonical_only = TRUE, 
                                                seen_canonical = seen_canonical,
                                                track_set = realizations_set)
    realization_count <- result$count
  } else {
    result <- stream_realizations_with_tracking(seed_chars, 1, character(length(seed_chars)), 
                                                conn_real, iupac_map, canonical_only = FALSE,
                                                track_set = realizations_set)
    realization_count <- result$count
  }
  
  close(conn_real)
  
  if (verbose) {
    cat("  Total realizations written:", realization_count, "\n\n")
  }
  
  # Step 2: Stream Hamming-1 neighbors
  if (verbose) {
    cat("Generating Hamming-1 neighbors to:", hamming_file, "\n")
    if (exclude_realizations_from_hamming) {
      cat("  (Excluding sequences that are already realizations)\n")
    }
  }
  
  if (deduplicate_hamming) {
    # Use hash table for deduplication (more memory but no duplicates)
    hamming_count <- generate_hamming1_deduplicated(realizations_file, hamming_file, 
                                                    seed_chars, non_n_positions, verbose,
                                                    canonical_only, realizations_set)
  } else {
    # Simple streaming (less memory but may have duplicates)
    counts <- process_realizations_chunked(seed_chars, non_n_positions, 
                                           realizations_file, hamming_file, 
                                           verbose = verbose, canonical_only = canonical_only,
                                           realizations_set = realizations_set)
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
    hamming1_count = hamming_count,
    canonical_only = canonical_only
  ))
}

# Function to generate deduplicated Hamming-1 neighbors using a hash table
generate_hamming1_deduplicated <- function(realizations_file, hamming_file, 
                                           seed_chars, non_n_positions, verbose = TRUE,
                                           canonical_only = FALSE, realizations_set = NULL) {
  bases <- c('A', 'C', 'G', 'T')
  seen_sequences <- new.env(hash = TRUE)  # Use environment as hash table
  
  conn_in <- file(realizations_file, "r")
  conn_out <- file(hamming_file, "w")
  
  total_written <- 0
  total_processed <- 0
  total_excluded <- 0
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
            
            # Get canonical form if requested
            if (canonical_only) {
              new_seq_str <- get_canonical(new_seq_str)
            }
            
            # Check if this sequence is a realization (and should be excluded)
            if (!is.null(realizations_set) && !is.null(realizations_set[[new_seq_str]])) {
              total_excluded <- total_excluded + 1
              next  # Skip this sequence
            }
            
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
            "unique Hamming-1 neighbors:", total_written)
        if (!is.null(realizations_set)) {
          cat(" (excluded", total_excluded, "that were realizations)")
        }
        cat("\n")
      }
    }
  }
  
  close(conn_in)
  close(conn_out)
  
  if (verbose && !is.null(realizations_set)) {
    cat("  Total sequences excluded (already realizations):", total_excluded, "\n")
  }
  
  return(total_written)
}

# Function to process multiple seeds with streaming
process_multiple_seeds_streaming <- function(seeds, output_dir = ".", verbose = TRUE,
                                             deduplicate_hamming = TRUE, canonical_only = FALSE,
                                             exclude_realizations_from_hamming = TRUE) {
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
    
    results[[i]] <- process_seed_streaming(seeds[i], output_prefix, verbose, 
                                           deduplicate_hamming, canonical_only,
                                           exclude_realizations_from_hamming)
  }
  
  names(results) <- seeds
  return(results)
}



# Example usage and main function
library("readr")
library("dplyr")

args <- commandArgs(trailingOnly = TRUE)
arrays=as.numeric(args[1]) #1-364
addition=as.numeric(args[2])
arrays=arrays+addition

#metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")
metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635

#sort based on the seed length

metadata=metadata[order(metadata$length),]

#length=10 
length=1
start_ind=(arrays-1)*length+1 #100
end_ind=arrays*length #100



if(end_ind > nrow(metadata)){
  end_ind=nrow(metadata)
}


path_local_scratch=Sys.getenv("LOCAL_SCRATCH")

lustre_folder="/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/"

for(i in start_ind:end_ind){

  print(i)
  motif= metadata$ID[i]
  seed=metadata$seed[i] 
  
  # Example with the provided seed
  

  #non_canon <- process_seed_streaming(seed1, output_prefix=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif), 
  #                                  verbose = TRUE, deduplicate_hamming = TRUE, canonical_only=FALSE)
  
  #theory <- theoretical_counts(seed, canonical = TRUE)
  #print(str(theory))
  
  canon <- process_seed_streaming(seed, output_prefix=paste0(path_local_scratch,"/", motif), 
                                    verbose = TRUE, deduplicate_hamming = TRUE, canonical_only=TRUE)
  print(str(canon))
  
  system(paste0("cp ", path_local_scratch,"/",motif,"_canonical_realizations.txt ", lustre_folder ))
  system(paste0("cp ", path_local_scratch,"/",motif,"_canonical_hamming1.txt ", lustre_folder ))

}

#comparison <- compare_theoretical_actual(seed)


#s=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(seed)
#motif=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(ID)


# Example usage and main function
# main <- function() {
#   # Example with the provided seed
#   cat("### EXAMPLE 1: Single seed with streaming ###\n")
#   seed1 <- "NNCCGGNNNNNNCCGGNN"
#   result1 <- process_seed_streaming(seed1, output_prefix = "example1", 
#                                     verbose = TRUE, deduplicate_hamming = TRUE,
#                                     canonical_only = FALSE)
#   
#   # Example with canonical k-mers only
#   cat("\n### EXAMPLE 2: Same seed with canonical k-mers only ###\n")
#   result2 <- process_seed_streaming(seed1, output_prefix = "example1_canonical", 
#                                     verbose = TRUE, deduplicate_hamming = TRUE,
#                                     canonical_only = TRUE)
#   
#   # Example with multiple seeds
#   cat("\n### EXAMPLE 3: Multiple seeds with canonical k-mers ###\n")
#   seeds <- c(
#     "NNCCGGNNNNNNCCGGNN",
#     "ATNGCN",
#     "RRYYWS"
#   )
#   
#   all_results <- process_multiple_seeds_streaming(
#     seeds, 
#     output_dir = "kmer_output",
#     verbose = TRUE,
#     deduplicate_hamming = TRUE,
#     canonical_only = TRUE
#   )
#   
#   # Summary
#   cat("\n### SUMMARY ###\n")
#   cat("----------------------------------------\n")
#   for (seed in names(all_results)) {
#     cat("Seed:", seed, "\n")
#     cat("  Canonical mode:", all_results[[seed]]$canonical_only, "\n")
#     cat("  Realizations file:", all_results[[seed]]$realizations_file, "\n")
#     cat("  Realizations count:", all_results[[seed]]$realization_count, "\n")
#     cat("  Hamming-1 file:", all_results[[seed]]$hamming1_file, "\n")
#     cat("  Hamming-1 count:", all_results[[seed]]$hamming1_count, "\n")
#     cat("----------------------------------------\n")
#   }
#   
#   return(all_results)
# }
# 
# # Function to process seeds from an input file
# process_seeds_from_file <- function(input_file, output_dir = "output", 
#                                    verbose = TRUE, deduplicate_hamming = TRUE,
#                                    canonical_only = FALSE) {
#   if (!file.exists(input_file)) {
#     stop("Input file does not exist:", input_file)
#   }
#   
#   # Read seeds from file
#   seeds <- readLines(input_file)
#   seeds <- seeds[seeds != ""]  # Remove empty lines
#   
#   cat("Found", length(seeds), "seeds in", input_file, "\n\n")
#   
#   # Process all seeds
#   results <- process_multiple_seeds_streaming(seeds, output_dir, verbose, 
#                                              deduplicate_hamming, canonical_only)
#   
#   # Write summary file
#   summary_file <- file.path(output_dir, "processing_summary.txt")
#   sink(summary_file)
#   cat("Processing Summary\n")
#   cat("==================\n")
#   cat("Date:", date(), "\n")
#   cat("Input file:", input_file, "\n")
#   cat("Canonical mode:", canonical_only, "\n")
#   cat("Total seeds processed:", length(seeds), "\n\n")
#   
#   for (seed in names(results)) {
#     cat("Seed:", seed, "\n")
#     cat("  Realizations:", results[[seed]]$realization_count, "\n")
#     cat("  Hamming-1 neighbors:", results[[seed]]$hamming1_count, "\n")
#     cat("\n")
#   }
#   sink()
#   
#   cat("\nSummary written to:", summary_file, "\n")
#   
#   return(results)
# }
# 
# # Run if script is executed directly
# if (!interactive()) {
#   # Check command line arguments
#   args <- commandArgs(trailingOnly = TRUE)
#   
#   if (length(args) > 0) {
#     # Process seeds from input file
#     input_file <- args[1]
#     output_dir <- ifelse(length(args) > 1, args[2], "output")
#     results <- process_seeds_from_file(input_file, output_dir)
#   } else {
#     # Run example
#     results <- main()
#   }
# }