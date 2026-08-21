library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("tidyverse")
library(R.utils)
library(stringr)
library(tidyr)




metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635

#which(metadata$seed=="NNNACGANNNNNNTCGTNNN") #992
#which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3")



#sort based on the seed length

metadata=metadata[order(metadata$length),]


results_path="/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/" #

#Filter out those which do not have either SELEX signal or background
metadata=metadata %>% filter(!(is.na(CSC_SELEX_filename) | is.na(CSC_SELEX_background_filename))) #3596

missing_signal=c()
missing_background=c()


for(index in 1:nrow(metadata)){ #0-9

  #index=1
  print(index)
  data=metadata[index,] 
  seed=data$seed
    
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", data$ID, ".txt")), silent=TRUE)
  
  
      
  kmer_signal_count_file=paste0("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/", data$ID, "_signal.txt")
  kmer_background_count_file=paste0("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/", data$ID, "_background.txt")
      
  kmers_signal_counts=try(read_delim(kmer_signal_count_file, col_names=FALSE, delim="\t"),silent=TRUE) #1126
  if(unique(class(kmers_signal_counts)=="try-error")){
    missing_signal=c(missing_signal, data$ID)
  }
  kmers_background_counts=try(read_delim(kmer_background_count_file, col_names=FALSE, delim="\t"), silent=TRUE) #1126
  if(unique(class(kmers_background_counts)=="try-error")){
    missing_background=c(missing_background, data$ID)
  }    
    
}

write_delim(as.data.frame(missing_signal), file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/count_kmers_with_Ns/rerun_missing.txt",
            col_names = FALSE)
  