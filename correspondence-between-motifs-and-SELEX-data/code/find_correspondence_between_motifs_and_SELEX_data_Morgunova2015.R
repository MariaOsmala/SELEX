library("readr")
library("dplyr")
library("stringr")
rm(list=ls())
metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")


data_path="/scratch/project_2013895/SELEX/data/"

#File names


#Morgunova2015: TFDP1.fastq.gz


metadata_Morgunova2015 = metadata %>% filter(study=="Morgunova2015")
study="Morgunova2015"

metadata_Morgunova2015$clone="DBD"

metadata_Morgunova2015$SELEX_filename="/scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp/E2F2.fastq.gz"

metadata_Morgunova2015$SELEX_background_filename=NA


metadata_Morgunova2015$unique_background=FALSE
metadata_Morgunova2015$cycle_background=NA

  




lookup=c(CSC_SELEX_filename="SELEX_filename",    
         CSC_SELEX_background_filename="SELEX_background_filename")


metadata_Morgunova2015 <- metadata_Morgunova2015 %>%
  rename(all_of(lookup) )

metadata_Morgunova2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp/", 
                                             "https://a3s.fi/Morgunova2015/", 
                                             metadata_Morgunova2015$CSC_SELEX_filename)


metadata_Morgunova2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp/", 
                                             "https://a3s.fi/Morgunova2015/", 
                                             metadata_Morgunova2015$CSC_SELEX_background_filename)


tmp = metadata_Morgunova2015 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,cycle_background, unique_background,seed,CSC_SELEX_filename, 
         CSC_SELEX_background_filename, Allas_SELEX_filename, 
         Allas_SELEX_background_filename)


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Morgunova2015.tsv",
            delim="\t")


