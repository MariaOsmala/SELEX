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

#Xie2025:NRL_FLI1_3_YYIII_TTTCTT40NTTTA.fastq.gz

metadata_Xie2025 = metadata %>% filter(study=="Xie2025")
study="Xie2025"

ENA_metadata=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Xie2025_from_ENA/filereport_read_run_PRJEB66722.tsv", delim="\t")
ENA_metadata$scientific_name %>% table(useNA="always") #Homo sapiens all? 
ENA_metadata$library_strategy %>% table(useNA="always") #
ENA_metadata$library_source %>% table(useNA="always") #
ENA_metadata$library_selection %>% table(useNA="always") #
ENA_metadata$study_accession %>% unique()
ENA_metadata$sample_alias %>% head()
ENA_metadata$sample_title %>% head()
splits=strsplit(ENA_metadata$sample_title, " ")
df <- as.data.frame(do.call(rbind, splits), stringsAsFactors = FALSE)

df$V1 %>% table()
df$V2 %>% table()
df$V3 %>% table()

names(df)=c("ENA_experiment", "sample","ENA_symbol")  

ENA_metadata=cbind(ENA_metadata, df[,c(1,3)])

splits=strsplit(ENA_metadata$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata$filename=do.call(rbind, splits)[,6]
ENA_metadata$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata$filename)

ENA_metadata=ENA_metadata %>% filter(ENA_experiment=="CAP-SELEX")

SELEX_data_filenames=dir(paste0(data_path,study, "/submitted_ftp/"))
splits <- strsplit(SELEX_data_filenames, "_")
lengths <- sapply(splits, length)
table(lengths) #all of equal length

df <- as.data.frame(do.call(rbind, splits), stringsAsFactors = FALSE)
df$filename=SELEX_data_filenames


colnames(df) <- c("TF1", "TF2", "cycle", "batch","file.ending", "filename")
df$ligand=do.call(rbind,strsplit(df$file.ending, "\\."))[,1]
df$file.ending=NULL

df$symbol=NA
df$symbol[which(df$TF1==df$TF2)]=df$TF1[which(df$TF1==df$TF2)]
df$symbol[which(df$TF1!=df$TF2)]=paste0(df$TF1[which(df$TF1!=df$TF2)],"_", df$TF2[which(df$TF1!=df$TF2)])

df=df[,c("TF1", "TF2", "symbol","cycle","batch","ligand","filename")]

ENA_metadata=ENA_metadata %>% left_join(df, by="filename")
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)


table(ENA_metadata$ENA_symbol==ENA_metadata$symbol)
#FALSE  TRUE 
#88 14817
#There is some issue with sample_alias and sample_title
ind=which(ENA_metadata$ENA_symbol!=ENA_metadata$symbol)
test=ENA_metadata[ind,names(ENA_metadata) %in% c("experiment_title", "experiment_alias","sample_alias", "sample_title","ENA_experiment", "ENA_symbol", "filename", "filename_without_ending","TF1","TF2",
                        "symbol", "cycle","batch",  "ligand", "csc_filename") ]


# Zero cycle background ---------------------------------------------------

ENA_metadata_background=read_delim("/projappl/project_2013895/SELEX/download-data/Data/ENA_metadata_input_libraries.tsv", delim="\t")
ENA_metadata_background$motif_derived_Xie2025="NO"
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)
ENA_metadata$motif_derived="NO"

metadata_Xie2025$cycle %>% table() #This is always 3


metadata_Xie2025$CSC_SELEX_filename=NA
metadata_Xie2025$CSC_SELEX_background_filename=NA
metadata_Xie2025$unique_background=FALSE
metadata_Xie2025$cycle_background=0

for(i in 1:nrow(metadata_Xie2025)){
  #i=1
  print(i)
  symbol=metadata_Xie2025$symbol[i]
  
  experiment=metadata_Xie2025$experiment[i]
  ligand=metadata_Xie2025$ligand[i]
  batch=metadata_Xie2025$batch[i]
  cycle=as.numeric(metadata_Xie2025$cycle[i])
  
  ENA_metadata %>% filter(symbol==.env$symbol)
  

 
  
  
  
  #Signal
  if(length(ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename))==1){
    metadata_Xie2025$CSC_SELEX_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol  & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename)
    ENA_metadata <- ENA_metadata %>%
      mutate(
        motif_derived = if_else(symbol==.env$symbol  & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle, "YES", motif_derived)
      )
  }
  

    if(length(ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename))==1){ 
      metadata_Xie2025$CSC_SELEX_background_filename[i]=ENA_metadata_background %>% filter(ligand==.env$ligand) %>% pull(csc_filename)
      
      ENA_metadata_background <- ENA_metadata_background %>%
        mutate(
          motif_derived = if_else(ligand==.env$ligand, "YES", motif_derived),
          motif_derived_Xie2025 = if_else(ligand==.env$ligand, "YES", motif_derived_Xie2025)
        )
    }
  
}




metadata_Xie2025$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/", 
                                             "https://a3s.fi/Xie2025/", 
                                             metadata_Xie2025$CSC_SELEX_filename)


metadata_Xie2025$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/", 
                                             "https://a3s.fi/Xie2025/", 
                                             metadata_Xie2025$CSC_SELEX_background_filename)

metadata_Xie2025$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/input_libraries/submitted_ftp/", 
                                                        "https://a3s.fi/input_libraries/", 
                                                        metadata_Xie2025$Allas_SELEX_background_filename)

metadata_Xie2025=metadata_Xie2025 %>% left_join( ENA_metadata %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_signal")), 
                                                     by=c("CSC_SELEX_filename"="csc_filename_signal"))

metadata_Xie2025=metadata_Xie2025 %>% left_join( ENA_metadata_background %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_background")), 
                                                     by=c("CSC_SELEX_background_filename"="csc_filename_background"))


metadata_Xie2025$study_accession_signal %>% unique()

missing_ind=which(is.na(metadata_Xie2025$CSC_SELEX_filename))
tmp=metadata_Xie2025[missing_ind,]

# VSX2\_TBX6\_TTTCGG40NAGA\_YAAII\_NAGGTGTTAATTN\_1\_3  
# ARNT\_TGGTCG40NTAA\_YT\_NGTCACGTGACN\_1\_3           
# ESR2\_TTCATA40NTCA\_YT\_NGGTCANNNNNNNNTGACCN\_1\_3  
# FOXK2\_TTGAAC40NTAGC\_YJI\_NRTAAAYAAAYAN\_1\_3       
# FOXM1\_TGACTG40NATG\_YT\_NGTAAACAAYRN\_1\_3         
# FOXM1\_TGACTG40NATG\_YT\_NRYAAAYAAACAN\_1\_3         
# FOXM1\_TGACTG40NATG\_YT\_NYGCATCAYAACRN\_1\_3        
# FOXP2\_TCGATT40NTTG\_YT\_NRCGTAAACAAN\_1\_3          
# HELT\_TAAGCC40NAGT\_YT\_NTCACGTGAC\_1\_3            
# OTP\_TAGATG40NTAT\_YPIII\_NTAATTRNNNNTAATTRN\_1\_3   
# SOX11\_TAGCCC40NGTA\_YT\_NACAATNNNNATTGTN\_1\_3     
# TCF23\_TGTTTA40NACG\_YT\_NGCCATTTGGTN\_1\_3          
# THHEX\_TCGAG40NCATT\_YT\_NCAATTNNNNNNNNNNAATTGN\_1\_3


missing_ind=which(is.na(metadata_Xie2025$CSC_SELEX_background_filename))

tmp=metadata_Xie2025[missing_ind,]
#tmp$ID

#"THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3"

metadata_Xie2025 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)

# "VSX2_TBX6_TTTCGG40NAGA_YAAII_NAGGTGTTAATTN_1_3"   "ARNT_TGGTCG40NTAA_YT_NGTCACGTGACN_1_3"            "ESR2_TTCATA40NTCA_YT_NGGTCANNNNNNNNTGACCN_1_3"   
# "FOXK2_TTGAAC40NTAGC_YJI_NRTAAAYAAAYAN_1_3"        "FOXM1_TGACTG40NATG_YT_NGTAAACAAYRN_1_3"           "FOXM1_TGACTG40NATG_YT_NRYAAAYAAACAN_1_3"         
# "FOXM1_TGACTG40NATG_YT_NYGCATCAYAACRN_1_3"         "FOXP2_TCGATT40NTTG_YT_NRCGTAAACAAN_1_3"           "HELT_TAAGCC40NAGT_YT_NTCACGTGAC_1_3"             
# "OTP_TAGATG40NTAT_YPIII_NTAATTRNNNNTAATTRN_1_3"    "SOX11_TAGCCC40NGTA_YT_NACAATNNNNATTGTN_1_3"       "TCF23_TGTTTA40NACG_YT_NGCCATTTGGTN_1_3"          
# "THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3"

#For which symbols, the motif was derived
symbols_with_data=ENA_metadata %>% filter(motif_derived=="YES") %>% pull(symbol) %>% unique()
#For which symbols, motif was not derived
symbols_without_data=ENA_metadata %>% filter(motif_derived=="NO") %>% pull(symbol) %>% unique()
symbols_without_data[!which(symbols_without_data %in% symbols_with_data)]


test=metadata_Xie2025 %>% filter(!is.na(CSC_SELEX_filename)) %>% filter(is.na(fastq_ftp_signal)) #9
metadata_Xie2025 %>% filter(!is.na(CSC_SELEX_background_filename)) %>% filter(is.na(fastq_ftp_background)) %>% nrow() #0

tmp = metadata_Xie2025 %>% 
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

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Xie2025.tsv",
            delim="\t")
#write ENA_background file
write_delim(ENA_metadata_background,"/projappl/project_2013895/SELEX/download-data/Data/ENA_metadata_input_libraries.tsv",
            delim="\t")

metadata_Xie2025$CSC_SELEX_filename %>% unique() %>% length() #1087

metadata_Xie2025 %>% nrow() #1348


exp_multip=names(which(metadata_Xie2025$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Xie2025 %>% filter(CSC_SELEX_filename %in% exp_multip) #175

tmp$CSC_SELEX_filename %>% table() %>% table()
# 2   3   4   5 
# 135  37  12   1 

which( (tmp$CSC_SELEX_filename %>% table() )==5)
