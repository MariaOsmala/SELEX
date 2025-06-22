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

#Jolma2015: NFATC4_HOXA1_3_AAE_TGAACG40NCAA.fastq.gz


metadata_Jolma2015 = metadata %>% filter(study=="Jolma2015")
study="Jolma2015"

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

colnames(df) <- c("TF1", "TF2", "cycle", "batch", "file.ending", "filename")
df$ligand=do.call(rbind,strsplit(df$file.ending, "\\."))[,1]
df$file.ending=NULL

df$symbol=NA
df$symbol[which(df$TF1==df$TF2)]=df$TF1[which(df$TF1==df$TF2)]
df$symbol[which(df$TF1!=df$TF2)]=paste0(df$TF1[which(df$TF1!=df$TF2)],"_", df$TF2[which(df$TF1!=df$TF2)])

df=df[,c("TF1", "TF2", "symbol","cycle","batch","ligand","filename")]

metadata_Jolma2015$SELEX_filename=NA
metadata_Jolma2015$SELEX_previous_cycle_filename=NA
metadata_Jolma2015$SELEX_ZeroCycle_filename=NA
  
for(i in 1:nrow(metadata_Jolma2015)){
  #i=1
  print(i)
  symbol=metadata_Jolma2015$symbol[i]
  
  experiment=metadata_Jolma2015$experiment[i]
  ligand=metadata_Jolma2015$ligand[i]
  batch=metadata_Jolma2015$batch[i]
  cycle=as.numeric(metadata_Jolma2015$cycle[i])
  
  df %>% filter(symbol==.env$symbol)
  
  metadata_Jolma2015$SELEX_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                              df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% 
                                                select(filename)
                                              )
  metadata_Jolma2015$SELEX_ZeroCycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol=="ZeroCycle" & ligand==.env$ligand) %>% select(filename))
  metadata_Jolma2015$SELEX_previous_cycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1) %>% select(filename))
}








lookup=c(SELEX_filename="CSC_SELEX_filename",               
         SELEX_previous_cycle_filename="CSC_SELEX_previous_cycle_filename",
         SELEX_ZeroCycle_filename="CSC_SELEX_ZeroCycle_filename")

metadata_Jolma2015 <- metadata_Jolma2015 %>%
  rename(all_of(lookup) )

metadata_Jolma2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2015/", 
                                             metadata_Jolma2015$CSC_SELEX_filename)


metadata_Jolma2015$Allas_SELEX_previous_cycle_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2015/", 
                                             metadata_Jolma2015$CSC_SELEX_previous_cycle_filename)
metadata_Jolma2015$Allas_SELEX_ZeroCycle_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                                      "https://a3s.fi/Jolma2015/", 
                                                      metadata_Jolma2015$CSC_SELEX_ZeroCycle_filename)


missing_ind=which(str_detect(metadata_Jolma2015$CSC_SELEX_filename, "/character"))
metadata_Jolma2015$CSC_SELEX_filename[missing_ind]=NA
metadata_Jolma2015$CSC_SELEX_previous_cycle_filename[missing_ind]=NA

metadata_Jolma2015$Allas_SELEX_filename[missing_ind]=NA
metadata_Jolma2015$Allas_SELEX_previous_cycle_filename[missing_ind]=NA

missing_ind=which(str_detect(metadata_Jolma2015$CSC_SELEX_ZeroCycle_filename, "/character"))
metadata_Jolma2015$CSC_SELEX_ZeroCycle_filename[missing_ind]=NA
metadata_Jolma2015$Allas_SELEX_ZeroCycle_filename[missing_ind]=NA

metadata_Jolma2015 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)

#"EHF_HT-SELEX_TCTTGA20NGTG_AG_NACCCGGAAGTA_2_3"
#"ELF3_HT-SELEX_TGACCT20NCCA_AG_NACCCGGAAGTAN_2_4"
#"ELF4_HT-SELEX_TGACTC20NTCA_AG_AACCCGGAAGTR_2_3" 

tmp = metadata_Jolma2015 %>% 
  select(ID,symbol,clone,Lambert2018_families,ligand, batch,cycle, seed,Allas_SELEX_filename, 
         Allas_SELEX_previous_cycle_filename, Allas_SELEX_ZeroCycle_filename)

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2015_Allas.tsv",
            delim="\t")


tmp = metadata_Jolma2015 %>% 
  select(ID,symbol,clone,Lambert2018_families,ligand, batch,cycle,seed,CSC_SELEX_filename, 
         CSC_SELEX_previous_cycle_filename, CSC_SELEX_ZeroCycle_filename)


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2015_CSC.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Jolma2015$CSC_SELEX_filename %>% unique() %>% length() #567

metadata_Jolma2015 %>% nrow() #820


exp_multip=names(which(metadata_Jolma2015$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Jolma2015 %>% filter(CSC_SELEX_filename %in% exp_multip) #175

tmp = tmp %>% filter(CSC_SELEX_filename=="/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/TBX1_TACGGA40NGGC_AI_4.fastq.gz")

tmp$CSC_SELEX_filename %>% table() %>% table()
# 2    3    4   5   6 
# 113  53   5   3   1 

which( (tmp$CSC_SELEX_filename %>% table() )==5)
