library("readr")
library("dplyr")
library("stringr")
rm(list=ls())

#SELEX data was downloaded from ENA in originally submitted format based on the following accession codes: 

# input libraries PRJEB20112
# Jolma 2013 PRJEB3289 study_accession
# Jolma 2015 PRJEB7934
# Morgunova 2015 PRJEB8671
# Nitta 2015 PRJEB7373
# Xie 2025 PRJEB66722
# Yin2017 PRJEB9797

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

strange_cycle=c("2b0","2u","3b0", "3b1", "3b1u",   "3u", "4b0",   "4b2",   "4u")

test=metadata %>% filter(cycle %in% strange_cycle) %>% select(ID, study, cycle)

test$study %>% table()
#Jolma2015 Nitta2015   Yin2017 
#83         2        55 

data_path="/scratch/project_2013895/SELEX/data/"

#File names

#Nitta2015: AR_TCGCAA40NTGT_1.fastq.gz
#Morgunova2015: TFDP1.fastq.gz

#Does this contain data for: "PAX3_HT-SELEX_TTAGGG20NGGA_AL_GTCACGCNNMATTAN_1_3", YES!

metadata_Nitta2015 = metadata %>% filter(study=="Nitta2015"| 
                                           ID=="PAX3_HT-SELEX_TTAGGG20NGGA_AL_GTCACGCNNMATTAN_1_3" | 
                                           ID=="NR1I2_HT-SELEX_TACTGG40NGGA_KR_RGTTCRNNNRGTTC_1_4" )
study="Nitta2015"

ENA_metadata=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Nitta2015_from_ENA/filereport_read_run_PRJEB7373.tsv")

ENA_metadata=ENA_metadata %>% filter(scientific_name=="Homo sapiens")
ENA_metadata$library_source %>% table(useNA="always") #
ENA_metadata$library_selection %>% table(useNA="always") #
ENA_metadata$study_alias %>% table(useNA="always") #
ENA_metadata$sample_alias
ENA_metadata$sample_title
  
ENA_metadata$submitted_ftp %>% head()

splits=strsplit(ENA_metadata$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata$filename=do.call(rbind, splits)[,6]
ENA_metadata$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata$filename)  

ENA_metadata %>% filter(sample_title=="PAX3") %>%pull(sample_alias)

ENA_metadata %>% filter(sample_title %in% c("FOS", "FOXA1","IRF2","NR1D2","NR1I2","RORB","RORB","SOX17")) %>% select(sample_title,sample_alias)   
#"NR1I2_HT-SELEX_TACTGG40NGGA_KR_RGTTCRNNNRGTTC_1_4" 

ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)
ENA_metadata$motif_derived="NO"

metadata_Nitta2015$CSC_SELEX_filename=NA
metadata_Nitta2015$CSC_SELEX_background_filename=NA

metadata_Nitta2015$unique_background=FALSE
metadata_Nitta2015$cycle_background=NA

  
for(i in 1:nrow(metadata_Nitta2015)){
  #i=1
  print(i)
  symbol=metadata_Nitta2015$symbol[i]
  
  ligand=metadata_Nitta2015$ligand[i]
  batch=metadata_Nitta2015$batch[i]
  cycle=metadata_Nitta2015$cycle[i]
  
  cycle_background=NULL
  unique_background=FALSE
  #is cycle on of the strange ones
  if(cycle %in% strange_cycle){ #2b0, 3b0,4b0 or 4u 
    
    if(length(grep("b", cycle))==1){
      cycles=strsplit(cycle, "b")[[1]]
      cycle=as.numeric(cycles[1]) #3u 2u 3b1u
      #Is there also u
      if(length(grep("u", cycles[2]))==1){
        cycles=strsplit(cycles[2], "u")[[1]]
        cycle_background=as.numeric(cycles)
        unique_background=TRUE
      }else{
        cycle_background=as.numeric(cycles[2])  
      }
      
      
    }else{
      #4u
      cycle=as.numeric(strsplit(cycle, "u")[[1]][1])
      cycle_background=cycle
      unique_background=TRUE
    }
    
  }else{
    cycle=as.numeric(cycle)
    cycle_background=cycle-1
  }
  
  sample_alias=paste0(batch, "_", ligand, "_", cycle )
  sample_alias_background=paste0(batch, "_", ligand, "_", cycle_background )
  
  metadata_Nitta2015$cycle_background[i]=cycle_background
  #filereport %>% filter(library_name==.env$symbol)
  
  
  if(length(ENA_metadata %>% filter(library_name==.env$symbol &  sample_alias==.env$sample_alias ) %>% pull(csc_filename))==1){
    metadata_Nitta2015$CSC_SELEX_filename[i]=ENA_metadata %>% filter(library_name==.env$symbol &  sample_alias==.env$sample_alias ) %>% pull(csc_filename)
  }
  if(length(ENA_metadata %>% filter(library_name==.env$symbol &  sample_alias==.env$sample_alias_background ) %>% pull(csc_filename) )==1){
    metadata_Nitta2015$CSC_SELEX_background_filename[i]=ENA_metadata %>% filter(library_name==.env$symbol &  sample_alias==.env$sample_alias_background ) %>% pull(csc_filename)
  }
  
  ENA_metadata <- ENA_metadata %>%
    mutate(
      motif_derived = if_else(library_name==.env$symbol &  sample_alias==.env$sample_alias , "YES", motif_derived)
    )
  
  ENA_metadata <- ENA_metadata %>%
    mutate(
      motif_derived = if_else(library_name==.env$symbol &  sample_alias==.env$sample_alias_background, "YES", motif_derived)
    )
  
  
  
}


metadata_Nitta2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp/", 
                                             "https://a3s.fi/Nitta2015/", 
                                             metadata_Nitta2015$CSC_SELEX_filename)


metadata_Nitta2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp/", 
                                             "https://a3s.fi/Nitta2015/", 
                                             metadata_Nitta2015$CSC_SELEX_background_filename)


#Add also other relevant ENA_metadata info to metadata_Jolm2013

metadata_Nitta2015=metadata_Nitta2015 %>% left_join( ENA_metadata %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_signal")), 
                                                     by=c("CSC_SELEX_filename"="csc_filename_signal"))

metadata_Nitta2015=metadata_Nitta2015 %>% left_join( ENA_metadata %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_background")), 
                                                     by=c("CSC_SELEX_background_filename"="csc_filename_background"))

#For which symbols, the motif was derived
symbols_with_data=ENA_metadata %>% filter(motif_derived=="YES") %>% pull(sample_title) %>% unique()
#For which symbols, motif was not derived
symbols_without_data=ENA_metadata %>% filter(motif_derived=="NO") %>% pull(sample_title) %>% unique()
symbols_without_data[!which(symbols_without_data %in% symbols_with_data)]

missing_ind=which(is.na(metadata_Nitta2015$CSC_SELEX_filename))

temp=metadata_Nitta2015[missing_ind,]
temp$ID

# SREBF2\_HT-SELEX\_TGAGAT20NGA\_Y\_NRTCACGCCAYN\_1\_3
# ONECUT2\_HT-SELEX\_TGGGCG30NCGT\_AH\_NNCGATCRATAWNN\_1\_2
# ONECUT1\_HT-SELEX\_TAGCTC20NTCT\_Y\_CRATCRATAWN\_1\_3  

missing_ind=which(is.na(metadata_Nitta2015$CSC_SELEX_background_filename))

#background also missing

metadata_Nitta2015 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)
metadata_Nitta2015 %>% filter(is.na(CSC_SELEX_background_filename)) %>% pull(ID)


tmp = metadata_Nitta2015 %>% 
  select(ID,motif_ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,cycle_background, unique_background,seed,CSC_SELEX_filename, 
         CSC_SELEX_background_filename, Allas_SELEX_filename, 
         Allas_SELEX_background_filename, 
         run_accession_signal,                 
         study_accession_signal,               secondary_study_accession_signal,      sample_accession_signal,              
         secondary_sample_accession_signal,     experiment_accession_signal,           submission_accession_signal,
         read_count_signal, base_count_signal,
         fastq_ftp_signal,                      submitted_ftp_signal,                  sra_ftp_signal,                       
         run_accession_background,              study_accession_background,            secondary_study_accession_background, 
         sample_accession_background,           secondary_sample_accession_background, experiment_accession_background,      
         submission_accession_background,       read_count_background, base_count_background,fastq_ftp_background,                  submitted_ftp_background,             
         sra_ftp_background                   )


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Nitta2015.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Nitta2015$CSC_SELEX_filename %>% unique() %>% length() #6

metadata_Nitta2015 %>% nrow() #10


exp_multip=names(which(metadata_Nitta2015$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Nitta2015 %>% filter(CSC_SELEX_filename %in% exp_multip) #


tmp$CSC_SELEX_filename %>% table() %>% table()

