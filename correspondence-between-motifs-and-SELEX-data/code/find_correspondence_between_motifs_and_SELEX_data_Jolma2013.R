library("readr")
library("dplyr")
library("stringr")

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

metadata_Jolma2013$SELEX_filename=NA
metadata_Jolma2013$SELEX_previous_cycle_filename=NA
metadata_Jolma2013$SELEX_ZeroCycle_filename=NA
  
for(i in 1:nrow(metadata_Jolma2013)){
  #i=1
  print(i)
  symbol=metadata_Jolma2013$symbol[i]
  
  experiment=metadata_Jolma2013$experiment[i]
  ligand=metadata_Jolma2013$ligand[i]
  batch=metadata_Jolma2013$batch[i]
  cycle=as.numeric(metadata_Jolma2013$cycle[i])
  
  df %>% filter(symbol==.env$symbol)
  
  metadata_Jolma2013$SELEX_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% select(filename))
  metadata_Jolma2013$SELEX_ZeroCycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol=="ZeroCycle" & ligand==.env$ligand) %>% select(filename))
  metadata_Jolma2013$SELEX_previous_cycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1) %>% select(filename))
}





lookup=c(CSC_SELEX_filename="SELEX_filename",            
         CSC_SELEX_previous_cycle_filename="SELEX_previous_cycle_filename",
         CSC_SELEX_ZeroCycle_filename="SELEX_ZeroCycle_filename")




metadata_Jolma2013 <- metadata_Jolma2013 %>%
  rename(all_of(lookup) )

metadata_Jolma2013$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2013/", 
                                             metadata_Jolma2013$CSC_SELEX_filename)


metadata_Jolma2013$Allas_SELEX_previous_cycle_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2013/", 
                                             metadata_Jolma2013$CSC_SELEX_previous_cycle_filename)
metadata_Jolma2013$Allas_SELEX_ZeroCycle_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/", 
                                                      "https://a3s.fi/Jolma2013/", 
                                                      metadata_Jolma2013$CSC_SELEX_ZeroCycle_filename)


missing_ind=which(str_detect(metadata_Jolma2013$CSC_SELEX_filename, "/character"))
metadata_Jolma2013$CSC_SELEX_filename[missing_ind]=NA
metadata_Jolma2013$CSC_SELEX_previous_cycle_filename[missing_ind]=NA

metadata_Jolma2013$Allas_SELEX_filename[missing_ind]=NA
metadata_Jolma2013$Allas_SELEX_previous_cycle_filename[missing_ind]=NA

missing_ind=which(str_detect(metadata_Jolma2013$CSC_SELEX_ZeroCycle_filename, "/character"))
metadata_Jolma2013$CSC_SELEX_ZeroCycle_filename[missing_ind]=NA
metadata_Jolma2013$Allas_SELEX_ZeroCycle_filename[missing_ind]=NA

metadata_Jolma2013 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)

#"EHF_HT-SELEX_TCTTGA20NGTG_AG_NACCCGGAAGTA_2_3"
#"ELF3_HT-SELEX_TGACCT20NCCA_AG_NACCCGGAAGTAN_2_4"
#"ELF4_HT-SELEX_TGACTC20NTCA_AG_AACCCGGAAGTR_2_3" 

tmp = metadata_Jolma2013 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle, seed,Allas_SELEX_filename, 
         Allas_SELEX_previous_cycle_filename)

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2013_Allas.tsv",
            delim="\t")


tmp = metadata_Jolma2013 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,seed,CSC_SELEX_filename, 
         CSC_SELEX_previous_cycle_filename)


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2013_CSC.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Jolma2013$CSC_SELEX_filename %>% unique() %>% length() #567

metadata_Jolma2013 %>% nrow() #820


exp_multip=names(which(metadata_Jolma2013$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Jolma2013 %>% filter(CSC_SELEX_filename %in% exp_multip) #175

tmp = tmp %>% filter(CSC_SELEX_filename=="/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/TBX1_TACGGA40NGGC_AI_4.fastq.gz")

tmp$CSC_SELEX_filename %>% table() %>% table()
# 2    3    4   5   6 
# 113  53   5   3   1 

which( (tmp$CSC_SELEX_filename %>% table() )==5)
