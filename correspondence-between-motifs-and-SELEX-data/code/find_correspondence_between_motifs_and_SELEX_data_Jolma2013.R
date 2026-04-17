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
#Jolma2013: NOTO_TGCGTT30NTGC_AI_1.fastq.gz
#Jolma2015: NFATC4_HOXA1_3_AAE_TGAACG40NCAA.fastq.gz
#Nitta2015: AR_TCGCAA40NTGT_1.fastq.gz
#Morgunova2015: TFDP1.fastq.gz
#Yin2015:
#YY2_FL_1_KV_TAGAGA40NCGC.fastq.gz
#YY2_FL_Bis_3_KAM_40NATATAATATTATTTATAATAAA.fastq.gz

#Xie2025:NRL_FLI1_3_YYIII_TTTCTT40NTTTA.fastq.gz

metadata_Jolma2013 = metadata %>% filter(study=="Jolma2013")
study="Jolma2013"

# metadata from ENA
#ERP001824(study accession PRJEB3289)
part1=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Jolma2013_from_ENA/filereport_read_run_ERP001824.tsv", delim="\t")
#ERP001826(PRJEB3291
part2=read_delim("/projappl/project_2013895/SELEX/download-data/experiments/Jolma2013_from_ENA/filereport_read_run_ERP001826.tsv", delim="\t")

ENA_metadata=rbind(part1, part2) #2756

ENA_metadata$scientific_name %>% table(useNA="always") #Homo sapiens all? 
ENA_metadata$library_strategy %>% table(useNA="always") #
#ChIP-Seq    OTHER     <NA> 
#  30     2726        0 
ENA_metadata$library_source %>% table(useNA="always") #
#GENOMIC SYNTHETIC      <NA> 
#  30      2726         0 
ENA_metadata$library_selection %>% table(useNA="always") #
# ChIP       other unspecified        <NA> 
#  30        2715          11           0 

ENA_metadata$study_alias %>% table(useNA="always") #
#HT-SELEX_ChIP-seq_confirmation HT-SELEX_human_TF-specificities                            <NA> 
#  30                            2726                               0 

ENA_metadata=ENA_metadata %>% filter(study_alias=="HT-SELEX_human_TF-specificities") #2726

ENA_metadata$study_accession %>% table()

ENA_metadata$sample_alias %>% head()
# "CTCF_full_AJ_TAGCGA20NGCT_3"
# "Dlx2_TCGCCA20NCCT_AA_3"
# "EGR1_TACTAT20NATC_AA_3"
# "FOXG1_TGGTTC20NGA_AA_2"
# "FOXG1_TGGTTC20NGA_AA_4"     
# "HSF1_TACCAC20NACC_AA_2"

ENA_metadata$sample_title %>% head()

table(ENA_metadata$sample_alias==ENA_metadata$sample_title) #all true

table(ENA_metadata$sample_alias==ENA_metadata$library_name) #all true

ENA_metadata$submitted_ftp %>% head()

splits=strsplit(ENA_metadata$submitted_ftp, "/")
sapply(splits, length) %>% table() #All 6
ENA_metadata$filename=do.call(rbind, splits)[,6]
ENA_metadata$filename_without_ending=gsub(".fastq.gz", "", ENA_metadata$filename)

#Path in csc
SELEX_data_filenames=dir(paste0(data_path,study, "/submitted_ftp/")) #2726

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

df[which(!is.na(df$V5)),]

# Step 3: Identify rows with and without "full" in column V2
is_full <- df$V2 == "full"

# --- Case A: Rows with "full" — need to swap V3 and V4
df[is_full, c("V3", "V4")] <- df[is_full, c("V4", "V3")]

# --- Case B: Rows without "full" — shift V2–V4 to the right, insert NA in V2
df[!is_full, c("V5", "V4", "V3")] <- df[!is_full, c("V4", "V3", "V2")]
df[!is_full, "V2"] <- NA

colnames(df) <- c("symbol", "clone", "ligand", "batch", "file.ending")
df$cycle=do.call(rbind,strsplit(df$file.ending, "\\."))[,1]
df$file.ending=NULL
df$filename=SELEX_data_filenames

#add df to ENA_metadata

ENA_metadata=ENA_metadata %>% left_join(df, by="filename")
ENA_metadata$csc_filename=paste0(data_path, study,"/submitted_ftp/", ENA_metadata$filename)
ENA_metadata$motif_derived="NO"

#Does this data contain all cycles

ENA_metadata %>% pull(cycle) %>% table()
#0   1   2   3   4   5   6   7 
#386 520 547 547 547 100  69  10 

metadata_Jolma2013$CSC_SELEX_filename=NA
metadata_Jolma2013$CSC_SELEX_background_filename=NA

metadata_Jolma2013$unique_background=FALSE
metadata_Jolma2013$cycle_background=NA

  
for(i in 1:nrow(metadata_Jolma2013)){
  #i=1
  print(i)
  symbol=metadata_Jolma2013$symbol[i]
  
  experiment=metadata_Jolma2013$experiment[i]
  ligand=metadata_Jolma2013$ligand[i]
  batch=metadata_Jolma2013$batch[i]
  cycle=as.numeric(metadata_Jolma2013$cycle[i])
  metadata_Jolma2013$cycle_background[i]=cycle-1
  ENA_metadata %>% filter(symbol==.env$symbol) %>% select(filename)
  
  if(length(ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename))==1){
    metadata_Jolma2013$CSC_SELEX_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% pull(csc_filename)
  }
  if(length(ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1) %>% pull(csc_filename))){
   metadata_Jolma2013$CSC_SELEX_background_filename[i]=ENA_metadata %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1) %>% pull(csc_filename)
  }
  ENA_metadata <- ENA_metadata %>%
    mutate(
      motif_derived = if_else(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle, "YES", motif_derived)
    )
  
  ENA_metadata <- ENA_metadata %>%
    mutate(
      motif_derived = if_else(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1, "YES", motif_derived)
    )
  
 
  
 
  
  
  
}


metadata_Jolma2013$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2013/", 
                                             metadata_Jolma2013$CSC_SELEX_filename)


metadata_Jolma2013$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2013/", 
                                             metadata_Jolma2013$CSC_SELEX_background_filename)


#Add also other relevant ENA_metadata info to metadata_Jolm2013

metadata_Jolma2013=metadata_Jolma2013 %>% left_join( ENA_metadata %>% 
                                         select(
                                           c("run_accession", "study_accession", "secondary_study_accession",
                                             "sample_accession", "secondary_sample_accession", "experiment_accession",
                                             "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                         ) %>%
                                         rename_with(~ paste0(.x, "_signal")), 
                                       by=c("CSC_SELEX_filename"="csc_filename_signal"))

metadata_Jolma2013=metadata_Jolma2013 %>% left_join( ENA_metadata %>% 
                                                       select(
                                                         c("run_accession", "study_accession", "secondary_study_accession",
                                                           "sample_accession", "secondary_sample_accession", "experiment_accession",
                                                           "submission_accession", "read_count", "base_count","fastq_ftp", "submitted_ftp", "sra_ftp", "csc_filename")
                                                       ) %>%
                                                       rename_with(~ paste0(.x, "_background")), 
                                                     by=c("CSC_SELEX_background_filename"="csc_filename_background"))



metadata_Jolma2013 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)
metadata_Jolma2013 %>% filter(is.na(CSC_SELEX_background_filename)) %>% pull(ID)
#"EHF_HT-SELEX_TCTTGA20NGTG_AG_NACCCGGAAGTA_2_3"
#"ELF3_HT-SELEX_TGACCT20NCCA_AG_NACCCGGAAGTAN_2_4"
#"ELF4_HT-SELEX_TGACTC20NTCA_AG_AACCCGGAAGTR_2_3" 

#For which symbols, the motif was derived
symbols_with_data=ENA_metadata %>% filter(motif_derived=="YES") %>% pull(symbol) %>% unique()
#For which symbols, motif was not derived
symbols_without_data=ENA_metadata %>% filter(motif_derived=="NO") %>% pull(symbol) %>% unique()
symbols_without_data[!which(symbols_without_data %in% symbols_with_data)]

tmp = metadata_Jolma2013 %>% 
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


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2013.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Jolma2013$CSC_SELEX_filename %>% unique() %>% length() #567

metadata_Jolma2013 %>% nrow() #820


exp_multip=names(which(metadata_Jolma2013$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Jolma2013 %>% filter(CSC_SELEX_filename %in% exp_multip) #175



tmp$CSC_SELEX_filename %>% table() %>% table()
# 2    3    4   5   6 
# 113  53   5   3   1 

which( (tmp$CSC_SELEX_filename %>% table() )==5)

tmp = tmp %>% filter(CSC_SELEX_filename=="/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/TBX1_TACGGA40NGGC_AI_4.fastq.gz")
