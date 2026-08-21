library("readr")
file="/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_19032026.tsv" #3596


metadata=read_delim(file, delim="\t")

is.na(metadata$CSC_SELEX_filename) %>% table()
is.na(metadata$CSC_SELEX_background_filename) %>% table()

(metadata$CSC_SELEX_filename=="") %>% table()
(metadata$CSC_SELEX_background_filename=="") %>% table()

#Filter out the corrupted files

metadata$fastq_ftp_signal %>% is.na() %>% table()
metadata$fastq_ftp_background %>% is.na() %>% table()

metadata=metadata %>% filter(!is.na(fastq_ftp_signal)) #3587


unique_combinations=unique(metadata[, c("CSC_SELEX_filename","CSC_SELEX_background_filename","ligand" )]) #2672

write_delim(unique_combinations, file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments/unique_signal_and_background.tsv", 
            delim="\t")
