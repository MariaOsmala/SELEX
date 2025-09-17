library(dplyr)
library(readr)


# Load metadata
# 3569
metadata <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins.tsv", 
                       delim = "\t")

# Output file path
# output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/combined_kmers_with_metadata_better.tsv"
output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/combined_kmers_with_metadata_better_maxscore_scaled.tsv"
#output_file <- "/scratch/project_2013895/SELEX/combined_kmers_protein_sequences/test.tsv"

# Write header only once
write_lines("ID\tprotein_sequence_1\tprotein_sequence_2\tprotein_sequence_3\tprotein_sequence_4\tkmer\tscore\trank", output_file)

# Get list of input files
#kmer_dir <- "/scratch/project_2013895/SELEX/streamed_kmers"
kmer_dir <- "/scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore"
kmer_files <- list.files(kmer_dir, pattern = "*.tsv", full.names = TRUE) #3755

#consider only those experiments that were successfull
#3731
successfull_exps=read.table("/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/successfull_exps.txt")

keep_ind=which( gsub(".tsv", "", do.call(rbind, strsplit(kmer_files, "/"))[,6]) %in% successfull_exps$V1)

kmer_files=kmer_files[keep_ind]

print(paste0("Succesfull kmer-files: ", length(kmer_files))) #3719
# Extract motif ID from filename
extract_id <- function(filename) {
  tools::file_path_sans_ext(basename(filename))
}

included_motifs=c()

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
  }else{
    included_motifs=c(included_motifs, motif_id)
  }
  
  # Read k-mer file in chunks
  con <- file(file, "r")
  while (length(line <- readLines(con, n = 10000)) > 0) {
    
    split_lines <- strsplit(line, "\\s+")  # or use "\\t" for tab-separated
    col_counts <- sapply(split_lines, length)
    
    if (any(col_counts != 3)) {
      warning(paste("Invalid format in file:", file, "- stopping."))
      #close(con)
      #stop("Aborted due to bad format.")  # Fully stops the script
      # Alternatively: break to just skip this file and move to the next one
      break
    }
    
    kmer_data <- read_table(paste(line, collapse = "\n"), col_names = c("kmer", "score", "rank"))
    
    repeated_meta <- meta_row[rep(1, nrow(kmer_data)), c("ID", 
                                                         "protein sequence 1", 
                                                         "protein sequence 2", 
                                                         "protein sequence 3", 
                                                         "protein sequence 4")]
    
    combined <- bind_cols(repeated_meta, kmer_data)
    
    # Rename columns for writing
    names(combined) <- c("ID", "protein_sequence_1", "protein_sequence_2", 
                         "protein_sequence_3", "protein_sequence_4", 
                         "kmer", "score", "rank")
    
    # Append to output file
    write_delim(combined, output_file, append = TRUE, col_names = FALSE, delim = "\t")
  }
  close(con)
}

saveRDS(included_motifs, "/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs.RDS")
#included_motifs=readRDS("/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs.RDS") #3360
# write.table(included_motifs, 
# "/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs.txt",
# quote=FALSE, 
# sep="\t",
# row.names=FALSE,
# col.names=FALSE)
