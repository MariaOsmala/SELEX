library("readr")
library("dplyr")
library("stringr")
rm(list=ls())
metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

strange_cycle=c("2b0","2u","3b0", "3b1", "3b1u",   "3u", "4b0",   "4b2",   "4u")

test=metadata %>% filter(cycle %in% strange_cycle) %>% select(ID, study, cycle)

test$study %>% table()
#Jolma2015 Nitta2015   Yin2017 
#83         2        55 

data_path="/scratch/project_2013895/SELEX/data/"

#SELEX data was downloaded from ENA in originally submitted format based on the following accession codes: 

# input libraries PRJEB20112
# Jolma 2013 PRJEB3289 study_accession
# Jolma 2015 PRJEB7934
# Morgunova 2015 PRJEB8671
# Nitta 2015 PRJEB7373
# Xie 2025 PRJEB66722
# Yin2017 PRJEB9797

#File names

#Jolma2015: NFATC4_HOXA1_3_AAE_TGAACG40NCAA.fastq.gz

metadata_Jolma2015 = metadata %>% filter(study=="Jolma2015" & ID!="PAX3_HT-SELEX_TTAGGG20NGGA_AL_GTCACGCNNMATTAN_1_3" &
                                          ID!= "NR1I2_HT-SELEX_TACTGG40NGGA_KR_RGTTCRNNNRGTTC_1_4" )
study="Jolma2015"

ENA_metadata=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Jolma2015_from_ENA/filereport_read_run_PRJEB7934.tsv", delim="\t")
ENA_metadata$scientific_name %>% table(useNA="always") #Homo sapiens all? 
ENA_metadata$library_strategy %>% table(useNA="always") #
#OTHER SELEX  <NA> 
#  2 14845     0 
ENA_metadata$library_source %>% table(useNA="always") #
#GENOMIC SYNTHETIC      <NA> 
#  2     14845         0 

ENA_metadata=ENA_metadata %>% filter(library_strategy=="SELEX")
ENA_metadata$library_selection %>% table(useNA="always") #
# other  <NA> 
# 14845     0 
ENA_metadata$study_alias %>% table(useNA="always") #
ENA_metadata$sample_alias %>% head()
ENA_metadata$sample_title %>% head()
splits=strsplit(ENA_metadata$sample_title, " ")
sapply(splits, length) %>% table() #Varies
lengths <- sapply(splits, length)
table(lengths)

max_len <- max(sapply(splits, length))

# Pad with NAs
padded <- lapply(splits, function(x) {
  length(x) <- max_len
  return(x)
})

df <- as.data.frame(do.call(rbind, padded), stringsAsFactors = FALSE)

df$V1 %>% table()
df$V2 %>% table()
df$V3 %>% table()

df[which(df$V1=="#N/A")]=#HT-SELEX# 

names(df)=c("ENA_experiment", "sample","ENA_symbol")  

ENA_metadata=cbind(ENA_metadata, df[,c(1,3)])

splits=strsplit(ENA_metadata$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata$filename=do.call(rbind, splits)[,6]
ENA_metadata$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata$filename)


SELEX_data_filenames=dir(paste0(data_path,study, "/submitted_ftp/"))
splits <- strsplit(SELEX_data_filenames, "_")
lengths <- sapply(splits, length)
table(lengths)

max_len <- max(sapply(splits, length))

# Pad with NAs
padded <- lapply(splits, function(x) {
  length(x) <- max_len
  return(x)
})

df <- as.data.frame(do.call(rbind, padded), stringsAsFactors = FALSE)
df$filename=SELEX_data_filenames
#Remove ChIP-seq?
df=df[-which(is.na(df$V5) ),]

#TF2 can be "htSELEX
colnames(df) <- c("TF1", "TF2", "cycle", "batch", "file.ending", "filename")
df$ligand=do.call(rbind,strsplit(df$file.ending, "\\."))[,1]
df$file.ending=NULL

df$symbol=NA
df$symbol[which(df$TF2=="htSELEX")]=df$TF1[which(df$TF2=="htSELEX")]

df$symbol[which(df$TF2!="htSELEX")]=paste0(df$TF1[which(df$TF2!="htSELEX")],"_", df$TF2[which(df$TF2!="htSELEX")])

df=df[,c("TF1", "TF2", "symbol","cycle","batch","ligand","filename")]

df$filename_experiment="CAP-SELEX"
df$filename_experiment[which(df$TF2=="htSELEX")]="HT-SELEX"

ENA_metadata=ENA_metadata %>% left_join(df, by="filename")
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)

table(ENA_metadata$ENA_experiment==ENA_metadata$filename_experiment)
#FALSE  TRUE 
#24 14821
ind=which(ENA_metadata$ENA_experiment!=ENA_metadata$filename_experiment)
test=ENA_metadata[ind,]
table(ENA_metadata$ENA_symbol==ENA_metadata$symbol)
#FALSE  TRUE 
#4 14817
ind=which(ENA_metadata$ENA_symbol!=ENA_metadata$symbol)
test=ENA_metadata[ind,]

# Zero cycle background ---------------------------------------------------

ENA_metadata_background=read_delim("/projappl/project_2013895/SELEX/download-data/Data/ENA_metadata_input_libraries.tsv", delim="\t")
ENA_metadata_background$motif_derived_Jolma2015="NO"
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)
ENA_metadata$motif_derived="NO"

iter_ind=which(!is.na(metadata_Jolma2015$cycle))

metadata_Jolma2015$CSC_SELEX_filename=NA
metadata_Jolma2015$CSC_SELEX_background_filename=NA
metadata_Jolma2015$unique_background=FALSE
metadata_Jolma2015$cycle_background=NA

  
for(i in iter_ind){
  #i=1
  print(i)
  symbol=metadata_Jolma2015$symbol[i]
  
  experiment=metadata_Jolma2015$experiment[i]
  ligand=metadata_Jolma2015$ligand[i]
  batch=metadata_Jolma2015$batch[i]
  cycle=metadata_Jolma2015$cycle[i]
  
  ENA_metadata %>% filter(symbol==.env$symbol)
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
  
  metadata_Jolma2015$unique_background[i]=unique_background
  metadata_Jolma2015$cycle_background[i]=cycle_background
  
  
  
  #Signal
  if(length(ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename))==1){
    metadata_Jolma2015$CSC_SELEX_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol  & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename)
    ENA_metadata <- ENA_metadata %>%
      mutate(
        motif_derived = if_else(symbol==.env$symbol  & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle, "YES", motif_derived)
      )
  }
  
  #Background
  if(cycle_background!=0){
    if(length(ENA_metadata %>% filter(symbol==.env$symbol &  ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background) %>% pull(csc_filename))==1){  
      metadata_Jolma2015$CSC_SELEX_background_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background) %>% pull(csc_filename)
      ENA_metadata <- ENA_metadata %>%
        mutate(
          motif_derived = if_else(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background, "YES", motif_derived))
    } 
  }else{
    if(length(ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename))==1){ 
      metadata_Jolma2015$CSC_SELEX_background_filename[i]=ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename)
      
      ENA_metadata_background <- ENA_metadata_background %>%
        mutate(
          motif_derived = if_else(ligand==.env$ligand, "YES", motif_derived),
          motif_derived_Jolma2015 = if_else(ligand==.env$ligand, "YES", motif_derived_Jolma2015)
        )
    }
  } #background
  
  
  
}


metadata_Jolma2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2015/", 
                                             metadata_Jolma2015$CSC_SELEX_filename)


metadata_Jolma2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                                      "https://a3s.fi/Jolma2015/", 
                                                      metadata_Jolma2015$CSC_SELEX_background_filename)


metadata_Jolma2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/input_libraries/submitted_ftp/", 
                                                      "https://a3s.fi/input_libraries/", 
                                                      metadata_Jolma2015$Allas_SELEX_background_filename)


#Add also other relevant ENA_metadata info to metadata_Jolm2013

setdiff(names(ENA_metadata), names(ENA_metadata_background))
setdiff( names(ENA_metadata_background), names(ENA_metadata))


ENA_all=dplyr::bind_rows(ENA_metadata, ENA_metadata_background)



metadata_Jolma2015=metadata_Jolma2015 %>% left_join( ENA_all %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_signal")), 
                                                     by=c("CSC_SELEX_filename"="csc_filename_signal"))

metadata_Jolma2015=metadata_Jolma2015 %>% left_join( ENA_all %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_background")), 
                                                     by=c("CSC_SELEX_background_filename"="csc_filename_background"))



missing_ind=which(is.na(metadata_Jolma2015$CSC_SELEX_filename))


tmp=metadata_Jolma2015[missing_ind,]

tmp$ID
# [1] "CEBPG_ATF4_CAP-SELEX_TGCGTC40NTTA_AAB_NNATGAYGCAAT_1_3b0" "BACH1_HT-SELEX_TTCCCC20NCCC_AL_ATGACTCAT_1_NA"            "FOS_HT-SELEX_TGAACT40NAAG_KR_NGATGACGTCATCR_2_4"         
# [4] "FOXA1_HT-SELEX_TTCTAA40NAAT_KN_TRNGTAAACA_1_3b1"          "IRF2_HT-SELEX_TTGCCC40NCTC_AAF_NAANCGAAASYR_1_3"          "NR1D2_HT-SELEX_TGAATT40NTAA_KR_TRGGTYASTAGGTCA_2_3"      
# [7] "RORB_HT-SELEX_TTCGGG40NGAG_KS_AANTAGGTCAGTAGGTCA_2_4"     "RORB_HT-SELEX_TTCGGG40NGAG_KS_AWNTAGGTCATGACCTANWT_2_4"   "SOX17_HT-SELEX_TATGCT40NACT_KO_ACCGAACAAT_1_4b2"

tmp$symbol
missing_ind=which(is.na(metadata_Jolma2015$CSC_SELEX_background_filename))

tmp=metadata_Jolma2015[missing_ind,]

tmp$ID #same as above

metadata_Jolma2015%>% filter(!is.na(CSC_SELEX_filename)) %>% filter(is.na(fastq_ftp_signal)) %>% nrow() #0
metadata_Jolma2015 %>% filter(!is.na(CSC_SELEX_background_filename)) %>% filter(is.na(fastq_ftp_background)) %>% nrow() #0

#For which symbols, the motif was derived
symbols_with_data=ENA_metadata %>% filter(motif_derived=="YES") %>% pull(symbol) %>% unique()
#For which symbols, motif was not derived
symbols_without_data=ENA_metadata %>% filter(motif_derived=="NO") %>% pull(symbol) %>% unique()
symbols_without_data[!which(symbols_without_data %in% symbols_with_data)]

tmp = metadata_Jolma2015 %>% 
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

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2015.tsv",
            delim="\t")

#write ENA_background file
write_delim(ENA_metadata_background,"/projappl/project_2013895/SELEX/download-data/Data/ENA_metadata_input_libraries.tsv",
            delim="\t")



#Are the SELEX data listed several times

metadata_Jolma2015$CSC_SELEX_filename %>% unique() %>% length() #345

metadata_Jolma2015 %>% nrow() #592


exp_multip=names(which(metadata_Jolma2015$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Jolma2015 %>% filter(CSC_SELEX_filename %in% exp_multip) #175

tmp$CSC_SELEX_filename %>% table() %>% table()
#2  3  4  5 
#87 37 15  7

which( (tmp$CSC_SELEX_filename %>% table() )==5)
