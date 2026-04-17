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

#File names
#Yin2017:
#YY2_FL_1_KV_TAGAGA40NCGC.fastq.gz
#YY2_FL_Bis_3_KAM_40NATATAATATTATTTATAATAAA.fastq.gz

metadata_Yin2017 = metadata %>% filter(study=="Yin2017" & experiment=="HT-SELEX") #864
study="Yin2017"

# metadata from ENA

#SELEX data was downloaded from ENA in originally submitted format based on the following accession codes: 

# input libraries PRJEB20112
# Jolma 2013 PRJEB3289 study_accession
# Jolma 2015 PRJEB7934
# Morgunova 2015 PRJEB8671
# Nitta 2015 PRJEB7373
# Xie 2025 PRJEB66722
# Yin2017 PRJEB9797

ENA_metadata=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Yin2017_from_ENA/filereport_read_run_PRJEB9797.tsv", delim="\t") #7443

ENA_metadata$scientific_name %>% table(useNA="always")
#Homo sapiens        Mus musculus synthetic construct                <NA> 
#  14                  49                7380                   0 

ENA_metadata$library_strategy %>% table(useNA="always") #
#ATAC-seq Bisulfite-Seq      ChIP-Seq         OTHER         SELEX          <NA> 
#  6            11            41             5          7380             0 
ENA_metadata$library_source %>% table(useNA="always") #
#GENOMIC SYNTHETIC      <NA> 
#  63      7380         0 
ENA_metadata$library_selection %>% table(useNA="always") #
#ChIP other  <NA> 
#  46  7397     0 

ENA_metadata$study_alias %>% table(useNA="always") #

ENA_metadata=ENA_metadata %>% filter(library_strategy=="SELEX") #7380

ENA_metadata$sample_alias %>% head()
# "ALX1_FL_2_KX_TCCTTG40NGGT"   "ALX3_eDBD_1_KO_TGCATG40NATG" "ALX3_eDBD_2_KO_TGCTCC40NCCG" "ALX3_FL_3_KV_TAGCTG40NCTA"   "ALX3_FL_4_KV_TTGCGA40NGTA"  
# "ALX4_eDBD_2_KW_TCTATT40NCAT"

ENA_metadata$sample_title %>% head()
#"mHT-SELEX sample ALX1_FL"   "mHT-SELEX sample ALX3_eDBD" "HT-SELEX sample ALX3_eDBD"  "HT-SELEX sample ALX3_FL"    "mHT-SELEX sample ALX3_FL"  
#"HT-SELEX sample ALX4_eDBD"

splits=strsplit(ENA_metadata$sample_title, " sample ")
sapply(splits, length) %>% table() #All 2
test=do.call(rbind, splits)
test[,1] %>% table(useNA="always")

# bHT-SELEX
# High salt concentration HT-SELEX
# High salt concentration mHT-SELEX
# HT-SELEX
# mHT-SELEX 

do.call(rbind, strsplit(test[,2], "_"))[,1] %>% table(useNA="always") #symbol
do.call(rbind, strsplit(test[,2], "_"))[,2] %>% table(useNA="always") #clones
#eDBD   FL <NA> 
#5024 2356    0 

ENA_metadata$ENA_experiment=test[,1]
ENA_metadata$ENA_clone=do.call(rbind, strsplit(test[,2], "_"))[,2]
ENA_metadata$ENA_symbol=do.call(rbind, strsplit(test[,2], "_"))[,1]


ENA_metadata$library_name %>% table() #unspecificed

table(ENA_metadata$sample_alias==ENA_metadata$sample_title) #all false

ENA_metadata$submitted_ftp %>% head()

splits=strsplit(ENA_metadata$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata$filename=do.call(rbind, splits)[,6]
ENA_metadata$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata$filename)


#path in csc
SELEX_data_filenames=dir(paste0(data_path,study, "/submitted_ftp/")) #7452
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

df$V1 %>% table()
df=df %>% filter(!(V1 %in% c("ATAC-seq", "Mouse", "Rabbit")))

df=df %>% filter(V2 %in% c("eDBD", "FL"))

df$V3 %>% table()

#Bis  Nor 

df=df %>% filter(!(V3 %in% c("Bis", "Nor")))  

df=df[,c(1:5,10)]
colnames(df) <- c("symbol", "clone", "cycle", "batch", "file.ending", "filename")

df$ligand=do.call(rbind,strsplit(df$file.ending, "\\."))[,1]

df$file.ending=NULL

df=df[, c("symbol", "clone", "cycle", "batch", "ligand", "filename")]

#df$sample_alias=gsub(".fastq.gz","",df$filename)

#add df to ENA_metadata

ENA_metadata=ENA_metadata %>% left_join(df, by="filename")

#We need the info whether the experiment is HT-SELEX or Methyl-SELEX 
ENA_metadata$ENA_experiment %>% table()
# "bHT-SELEX"
# "High salt concentration HT-SELEX"
# "High salt concentration mHT-SELEX"
# "HT-SELEX" 
# "mHT-SELEX"


# Zero cycle background ---------------------------------------------------

ENA_metadata_background=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/input_libraries/filereport_read_run_ERP022234.tsv", delim="\t")

splits=strsplit(ENA_metadata_background$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata_background$filename=do.call(rbind, splits)[,6]
ENA_metadata_background$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata_background$filename)


input_library_filenames=dir(paste0(data_path,"input_libraries", "/submitted_ftp/")) #1029

df_input <- as.data.frame(do.call(rbind, strsplit(input_library_filenames, "_")))

df_input$filename=input_library_filenames
df_input$V1=NULL
df_input$V3=NULL
df_input$V4=NULL
names(df_input)[1]="ligand"

ENA_metadata_background=ENA_metadata_background %>% left_join(df_input, by="filename")

# Find corresponding filenames --------------------------------------------
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)
ENA_metadata_background$csc_filename=paste0(data_path, "input_libraries","/submitted_ftp/", ENA_metadata_background$filename)
ENA_metadata$motif_derived="NO"
ENA_metadata_background$motif_derived="NO"
metadata_Yin2017$CSC_SELEX_filename=NA
metadata_Yin2017$CSC_SELEX_background_filename=NA


metadata_Yin2017$unique_background=FALSE
metadata_Yin2017$cycle_background=NA

  
for(i in 1:nrow(metadata_Yin2017)){
  #i=1
  print(i)
  symbol=metadata_Yin2017$symbol[i]
  
  experiment=metadata_Yin2017$experiment[i] #HT-SELEX Methyl-HT-SELEX
  ligand=metadata_Yin2017$ligand[i]
  batch=metadata_Yin2017$batch[i]
  cycle=metadata_Yin2017$cycle[i]
  cycle_background=NULL
  unique_background=FALSE
  #is cycle on of the strange ones
  if(cycle %in% strange_cycle){ #2b0, 3b0,4b0 or 4u 
    
    if(length(grep("b", cycle))==1){
      cycle_background=0
      cycle=as.numeric(strsplit(cycle, "b")[[1]][1])
      
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
  
  metadata_Yin2017$unique_background[i]=unique_background
  metadata_Yin2017$cycle_background[i]=cycle_background
  
  tmp=ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment)
  
  
  #Signal
  if(length(ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename))==1){
    metadata_Yin2017$CSC_SELEX_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename)
    ENA_metadata <- ENA_metadata %>%
      mutate(
        motif_derived = if_else(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle, "YES", motif_derived)
      )
  }
  
  #Background
  if(cycle_background!=0){
    if(length(ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background) %>% pull(csc_filename))==1){  
      metadata_Yin2017$CSC_SELEX_background_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background) %>% pull(csc_filename)
      ENA_metadata <- ENA_metadata %>%
        mutate(
          motif_derived = if_else(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background, "YES", motif_derived)
        )
    } 
  }else{
    if(length(ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename))==1){ 
      metadata_Yin2017$CSC_SELEX_background_filename[i]=ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename)
      
      ENA_metadata_background <- ENA_metadata_background %>%
        mutate(
          motif_derived = if_else(ligand==.env$ligand, "YES", motif_derived)
        )
    }
  } #background
}

metadata_Yin2017$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp/", 
                                             "https://a3s.fi/Yin2017/", 
                                             metadata_Yin2017$CSC_SELEX_filename)


metadata_Yin2017$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp/", 
                                             "https://a3s.fi/Yin2017/", 
                                             metadata_Yin2017$CSC_SELEX_background_filename)


metadata_Yin2017$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/input_libraries/submitted_ftp/", 
                                                      "https://a3s.fi/input_libraries/", 
                                                      metadata_Yin2017$Allas_SELEX_background_filename)


missing_ind=which(is.na(metadata_Yin2017$CSC_SELEX_filename))
temp=metadata_Yin2017[missing_ind,]
#temp$symbol "ZSCAN5A could be ZSCAN5" "ZSCAN31 could be ZNF323"

temp$ID

missing_ind=which(is.na(metadata_Yin2017$CSC_SELEX_background_filename))
temp=metadata_Yin2017[missing_ind,]

metadata_Yin2017 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)

#"ZSCAN5A_HT-SELEX_TCGCCC40NCAT_KR_NYGTCCCYCCCCAAANMN_2_2"
#"ZSCAN31_HT-SELEX_TGGAGA40NCCA_KV_GCATAACKGCCCTGCKKCN_2_4" 

metadata_Yin2017 %>% filter(is.na(CSC_SELEX_background_filename)) %>% pull(ID)

#Add also other relevant ENA_metadata info to metadata_Jolm2013

setdiff(names(ENA_metadata), names(ENA_metadata_background))
setdiff( names(ENA_metadata_background), names(ENA_metadata))

ENA_metadata$fastq_bytes=as.numeric(ENA_metadata$fastq_bytes)
ENA_metadata_background$fastq_bytes=as.numeric(ENA_metadata_background$fastq_bytes)

ENA_metadata$submitted_bytes=as.numeric(ENA_metadata$submitted_bytes)
ENA_metadata_background$submitted_bytes=as.numeric(ENA_metadata_background$submitted_bytes)

ENA_all=dplyr::bind_rows(ENA_metadata, ENA_metadata_background)

metadata_Yin2017=metadata_Yin2017 %>% left_join( ENA_all %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession", "ENA_experiment",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_signal")), 
                                                     by=c("CSC_SELEX_filename"="csc_filename_signal"))

metadata_Yin2017=metadata_Yin2017 %>% left_join( ENA_all %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession","ENA_experiment",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_background")), 
                                                     by=c("CSC_SELEX_background_filename"="csc_filename_background"))



metadata_Yin2017%>% filter(!is.na(CSC_SELEX_filename)) %>% filter(is.na(fastq_ftp_signal)) %>% nrow() #0
metadata_Yin2017 %>% filter(!is.na(CSC_SELEX_background_filename)) %>% filter(is.na(fastq_ftp_background)) %>% nrow() #0



tmp = metadata_Yin2017 %>% 
  select(ID,motif_ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,cycle_background, unique_background,seed,CSC_SELEX_filename, 
         CSC_SELEX_background_filename, Allas_SELEX_filename, 
         Allas_SELEX_background_filename, 
         run_accession_signal,                 
         study_accession_signal,               secondary_study_accession_signal,      sample_accession_signal,              
         secondary_sample_accession_signal,     experiment_accession_signal,     ENA_experiment_signal,      submission_accession_signal,
         read_count_signal, base_count_signal,
         fastq_ftp_signal,                      submitted_ftp_signal,                  sra_ftp_signal,                       
         run_accession_background,              study_accession_background,            secondary_study_accession_background, 
         sample_accession_background,           secondary_sample_accession_background, experiment_accession_background,   ENA_experiment_background,   
         submission_accession_background,       read_count_background, base_count_background,fastq_ftp_background,                  submitted_ftp_background,             
         sra_ftp_background                   )




write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Yin2017.tsv",
            delim="\t")

#write ENA_background file
write_delim(ENA_metadata_background,"/projappl/project_2013895/SELEX/download-data/Data/ENA_metadata_input_libraries.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Yin2017$CSC_SELEX_filename %>% unique() %>% length() #683

metadata_Yin2017 %>% nrow() #864


exp_multip=names(which(metadata_Yin2017$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Yin2017 %>% filter(CSC_SELEX_filename %in% exp_multip) #358



tmp$CSC_SELEX_filename %>% table() %>% table()
#2   3   4 
#120  18   8

which( (tmp$CSC_SELEX_filename %>% table() )==4)
