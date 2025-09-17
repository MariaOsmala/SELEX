
library("readr")
library("tidyverse")
library("dplyr")
SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv")

SELEX_data %>% filter(ID=="MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2") %>% pull(CSC_SELEX_filename, CSC_SELEX_background_filename)


kmers=read_delim("/scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore/MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2.tsv",
                 col_names = FALSE)

as.character(kmers[1,1])
