library("readr")
library("dplyr")
library("tidyverse")
metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t")


metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))



missing=metadata %>% filter(is.na(CSC_SELEX_filename)| is.na(CSC_SELEX_background_filename)) %>% select(ID)

metadata=metadata %>% filter(!(ID %in% missing$ID))

#is.na(metadata$CSC_SELEX_filename) %>% table()
#is.na(metadata$CSC_SELEX_background_filename) %>% table()

#sort based on the seed length

metadata=metadata[order(metadata$length),]

#which(metadata$seed=="NNNACGANNNNNNTCGTNNN") #3180

metadata %>% filter(seed=="NNNACGANNNNNNTCGTNNN") %>% pull(ID)

write_delim(metadata, "/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv", delim="\t") #3573



