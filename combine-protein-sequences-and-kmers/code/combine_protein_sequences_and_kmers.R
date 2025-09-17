library(dplyr)
library(readr)


metadata <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins.tsv", 
                      delim="\t")

names(metadata)


process_kmer_file <- function(kmer_filepath, motif_id, metadata_df) {
  # Read the k-mer file (assuming whitespace-separated: kmer, score, rank)
  kmer_data <- read_table(kmer_filepath, col_names = c("kmer", "score", "rank"))
  
  # Get metadata for this motif
  meta_row <- metadata_df %>% filter(ID == motif_id)
  
  if (nrow(meta_row) == 0) {
    warning(paste("No metadata found for", motif_id))
    return(NULL)
  }
  
  # Repeat metadata across all k-mer rows
  repeated_meta <- meta_row[rep(1, nrow(kmer_data)), c("ID", 
                                                       "protein sequence 1", 
                                                       "protein sequence 2", 
                                                       "protein sequence 3", 
                                                       "protein sequence 4")]
  
  # Combine metadata and k-mer info
  combined <- bind_cols(repeated_meta, kmer_data)
  return(combined)
}


# Directory where k-mer files are stored
kmer_dir <- "/scratch/project_2013895/SELEX/streamed_kmers/"

# List all k-mer files
kmer_files <- list.files(kmer_dir, pattern = "*.tsv", full.names = TRUE)

# Extract IDs from filenames (adjust regex as needed)
extract_id <- function(filename) {
  # Assumes filenames like kmers/BCL6B_HT-SELEX_...tsv
  tools::file_path_sans_ext(basename(filename))
}

# Process all files and combine
all_kmers <- lapply(kmer_files, function(file) {
  motif_id <- extract_id(file)
  process_kmer_file(file, motif_id, metadata)
})

save.image("/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/all_kmers.RData")

# Combine into one dataframe
final_table <- bind_rows(all_kmers)

save.image("/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/final_table.RData")

write_delim(final_table, "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/combined_kmers_with_metadata.csv", row.names = FALSE, 
            delim="\t")

