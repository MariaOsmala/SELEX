library(dplyr)
library(readr)

rm(list=ls())
# Load metadata
# 3569
#metadata <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins.tsv", 
#                       delim = "\t")

#This is needed for motif IDs
metadata_motif_IDs <- read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv", delim="\t") 

metadata_motif_IDs <- metadata_motif_IDs %>% dplyr::select(ID, motif_ID)

#3594, lambda has been already added
metadata <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_04032026.tsv", 
                       delim = "\t") #3597

metadata=metadata %>%left_join(metadata_motif_IDs, by="ID") %>% relocate(motif_ID, .after=ID)


is.na(metadata$lambda) %>% table()
#206 TRUE, 3388 FALSE

(metadata$CSC_SELEX_filename=="") %>% table()
(metadata$CSC_SELEX_background_filename=="") %>% table()

is.na(metadata$CSC_SELEX_filename=="") %>% table()
is.na(metadata$CSC_SELEX_background_filename=="") %>% table()


# Output file path
# output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/combined_kmers_with_metadata_better.tsv"
#output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/combined_kmers_with_metadata_better_maxscore_scaled.tsv"
#output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/test.tsv"
output_file <- "/scratch/project_2013895/SELEX/sequence-to-affinity-data/all.tsv"
separate_folder <- "/scratch/project_2013895/SELEX/sequence-to-affinity-data/separate_with_protein_sequences/"
separate_folder_without_proteins <- "/scratch/project_2013895/SELEX/sequence-to-affinity-data/separate_without_protein_sequences/"

# Write header only once
write_lines("motif_ID\tprotein_sequence_1\tprotein_sequence_2\tprotein_sequence_3\tprotein_sequence_4\tkmer\t count_score\tcount_score_rank\tmotif_match_score\tmotif_match_score_rank", output_file)

# Get list of input files
#kmer_dir <- "/scratch/project_2013895/SELEX/streamed_kmers"
#kmer_dir <- "/scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore"
#kmer_files <- list.files(kmer_dir, pattern = "*.tsv", full.names = TRUE) #3755

kmer_dir <- "/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/"
kmer_files <- list.files(kmer_dir, pattern = "*.tsv", full.names = TRUE) #3469

# Extract motif ID from filename
extract_id <- function(filename) {
  tools::file_path_sans_ext(basename(filename))
}

included_motifs=c()

#Remove cases for which we do not have correct protein info
motifs_with_protein_info <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_protein_info_04032026.tsv", delim="\t") #3224




# Process files one-by-one
for (file in kmer_files) {
  #file=kmer_files[1]
  #print(strsplit(file, "/")[[1]][6])
  print(file)
  motif_id <- extract_id(file)
  
  meta_row <- metadata %>% filter(ID == motif_id)
  
  if (nrow(meta_row) == 0) {
    warning(paste("No metadata found for", motif_id))
    next
  }
  
  kmer_data=read_delim(file=file, delim="\t")
  if( names(table(is.na(kmer_data$Corrected_count_scaled)))=="TRUE"){
    next
  }
  
  if(!(motif_id %in% motifs_with_protein_info$ID)){
    #There is no protein info for this motif
    next
  }
  
  #Do not handle this case if TF1_copies is NA or TF1_copies or TF2_copies contain ?. This is probably handled already above

  if(is.na(meta_row$TF1_copies) | (grepl("\\?", meta_row$TF1_copies) | grepl("\\?", meta_row$TF2_copies)) ){
    next
  }
  
  repeated_meta <- meta_row[rep(1, nrow(kmer_data)), c("motif_ID",
                                                       "protein sequence 1", 
                                                       "protein sequence 2", 
                                                       "protein sequence 3", 
                                                       "protein sequence 4")]
  
  combined <- bind_cols(repeated_meta, kmer_data)
  
  # Rename columns for writing
  names(combined) <- c("motif_ID", "protein_sequence_1", "protein_sequence_2", 
                       "protein_sequence_3", "protein_sequence_4", 
                       "kmer", "count_score", "count_score_rank","motif_match_score","motif_match_score_rank")
  
  
  
 
  write_delim(combined, output_file, append = TRUE, col_names = FALSE, delim = "\t")
  #write also separately without ID
  write_delim(combined %>% dplyr::select(-motif_ID), paste0(separate_folder, meta_row$motif_ID, ".tsv"), append = FALSE, col_names = TRUE, delim = "\t")
  #write also separately without protein sequences
  write_delim(combined %>% dplyr::select(-motif_ID, -protein_sequence_1, -protein_sequence_2,-protein_sequence_3,-protein_sequence_4 ), paste0(separate_folder_without_proteins, meta_row$motif_ID, ".tsv"), append = FALSE, col_names = TRUE, delim = "\t")
              
 
  included_motifs=c(included_motifs, motif_id)
 
}

saveRDS(included_motifs, "/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs_final.RDS")

#included_motifs=readRDS("/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs_final.RDS") #3107
write.table(included_motifs, 
 "/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs_final.txt",
 quote=FALSE, 
 sep="\t",
 row.names=FALSE,
 col.names=FALSE)


#Write included motifs to metadata table

#table(metadata$ID %in% included_motifs)

#FALSE  TRUE 
#503  3091 

metadata$included_in_S2A=FALSE
metadata$included_in_S2A[metadata$ID %in% included_motifs]=TRUE
write_delim(metadata, "/projappl/project_2013895/motif_metadata/motifs_with_S2A_and_proteins.tsv", 
                                              delim = "\t")

# Read k-mer file in chunks
# con <- file(file, "r")
# while (length(line <- readLines(con, n = 10000)) > 0) {
#   
#   split_lines <- strsplit(line, "\\s+")  # or use "\\t" for tab-separated
#   col_counts <- sapply(split_lines, length)
#   
#   if (any(col_counts != 3)) {
#     warning(paste("Invalid format in file:", file, "- stopping."))
#     #close(con)
#     #stop("Aborted due to bad format.")  # Fully stops the script
#     # Alternatively: break to just skip this file and move to the next one
#     break
#   }
#   
#   kmer_data <- read_table(paste(line, collapse = "\n"), col_names = c("kmer", "score", "rank"))
#   
#   repeated_meta <- meta_row[rep(1, nrow(kmer_data)), c("ID", 
#                                                        "protein sequence 1", 
#                                                        "protein sequence 2", 
#                                                        "protein sequence 3", 
#                                                        "protein sequence 4")]
#   
#   combined <- bind_cols(repeated_meta, kmer_data)
#   
#   # Rename columns for writing
#   names(combined) <- c("ID", "protein_sequence_1", "protein_sequence_2", 
#                        "protein_sequence_3", "protein_sequence_4", 
#                        "kmer", "score", "rank")
#   
#   # Append to output file
#   write_delim(combined, output_file, append = TRUE, col_names = FALSE, delim = "\t")
# }
# close(con)

