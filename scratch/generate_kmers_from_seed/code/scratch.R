# start_ind=arrays*10+1 #100
# end_ind=(arrays+1)*10 #100
# len=10 #100
# 
# #3933 motifs
# if(end_ind>393){ # 398 
#   end_ind=3933
# }

# This code could also score the k-mers based on PWM
# sort the k-mers based on the scores and give rank 

# Function to expand IUPAC ----------------------------------------------

#iupac_map <- Biostrings::IUPAC_CODE_MAP


# stream_iupac_kmers_to_file <- function(iupac_kmer, output_file) {
#   con <- file(output_file, open = "w")
#   on.exit(close(con))
#   
#   iupac_map <- Biostrings::IUPAC_CODE_MAP
#   chars <- unlist(strsplit(iupac_kmer, ""))
#   
#   # Convert to list of vectors of valid bases
#   base_options <- lapply(chars, function(ch) strsplit(iupac_map[[ch]], "")[[1]])
#   
#   # Recursive streaming function
#   stream_recursive <- function(pos = 1, prefix = "") {
#     if (pos > length(base_options)) {
#       writeLines(prefix, con)
#     } else {
#       for (base in base_options[[pos]]) {
#         stream_recursive(pos + 1, paste0(prefix, base))
#       }
#     }
#   }
#   
#   stream_recursive()
# }

# Example usage: Expand "ATNGCR" into all valid k-mers and write to a file
# stream_iupac_kmers_to_file(metadata$seed[arrays], paste0("/scratch/project_2013895/SELEX/streamed_kmers/",metadata$ID[arrays],".txt"))

# count_possible_sequences <- function(iupac_seq) {
#   iupac_map <- Biostrings::IUPAC_CODE_MAP
#   chars <- unlist(strsplit(iupac_seq, ""))
#   
#   # Number of options per character
#   possibilities_per_position <- sapply(chars, function(ch) {
#     nchar(iupac_map[[ch]])
#   })
#   
#   # Multiply across all positions
#   total_possibilities <- prod(possibilities_per_position)
#   return(total_possibilities)
# }

# Example:
# count_possible_sequences(metadata$seed[arrays])  # Returns: 16777216

# stream_unique_canonical_kmers_to_file <- function(iupac_seq, output_file) {
#   iupac_map <- Biostrings::IUPAC_CODE_MAP
#   chars <- unlist(strsplit(iupac_seq, ""))
#   base_options <- lapply(chars, function(ch) strsplit(iupac_map[[ch]], "")[[1]])
#   
#   seen <- new.env(hash = TRUE, parent = emptyenv())
#   con <- file(output_file, open = "w")
#   on.exit(close(con))
#   
#   stream_recursive <- function(pos = 1, prefix = "") {
#     if (pos > length(base_options)) {
#       kmer <- DNAString(prefix)
#       revcomp <- as.character(reverseComplement(kmer))
#       canonical <- min(as.character(kmer), revcomp)
#       
#       if (!exists(canonical, envir = seen)) {
#         assign(canonical, TRUE, envir = seen)
#         writeLines(canonical, con)
#       }
#     } else {
#       for (base in base_options[[pos]]) {
#         stream_recursive(pos + 1, paste0(prefix, base))
#       }
#     }
#   }
#   
#   stream_recursive()
# }

#stream_unique_canonical_kmers_to_file(metadata$seed[arrays], paste0("/scratch/project_2013895/SELEX/streamed_kmers/",metadata$ID[arrays],"_canonical.txt"))
