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
#Yin2017:
#YY2_FL_1_KV_TAGAGA40NCGC.fastq.gz
#YY2_FL_Bis_3_KAM_40NATATAATATTATTTATAATAAA.fastq.gz

metadata_Yin2017 = metadata %>% filter(study=="Yin2017")
study="Yin2017"

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

df$sample_alias=gsub(".fastq.gz","",df$filename)

#We need the info whether the experiment is HT-SELEX or Methyl-SELEX 

SELEX_metadata=read_delim("/scratch/project_2013895/SELEX/data/Yin2017/filereport_read_run_PRJEB9797_tsv.txt")

df=df %>%
  left_join(SELEX_metadata %>% select(sample_alias, sample_title), by="sample_alias") 

df=cbind(df,do.call(rbind,strsplit(df$sample_title, " sample ")))
names(df)[(ncol(df)-1):ncol(df)]=c("experiment", "symbol_clone")

df$experiment=gsub("mHT-SELEX","Methyl-HT-SELEX", df$experiment)


# Zero cycle background ---------------------------------------------------

input_library_filenames=dir(paste0(data_path,"input_libraries", "/submitted_ftp/"))

df_input <- as.data.frame(do.call(rbind, strsplit(input_library_filenames, "_")))

df_input$filename=input_library_filenames
df_input$V1=NULL
df_input$V3=NULL
df_input$V4=NULL
names(df_input)[1]="ligand"

# Find corresponding filenames --------------------------------------------


metadata_Yin2017$SELEX_filename=NA
metadata_Yin2017$SELEX_previous_cycle_filename=NA
metadata_Yin2017$SELEX_ZeroCycle_filename=NA
  
for(i in 1:nrow(metadata_Yin2017)){
  #i=1
  print(i)
  symbol=metadata_Yin2017$symbol[i]
  
  experiment=metadata_Yin2017$experiment[i] #HT-SELEX Methyl-HT-SELEX
  ligand=metadata_Yin2017$ligand[i]
  batch=metadata_Yin2017$batch[i]
  cycle=metadata_Yin2017$cycle[i]
  zero_background=FALSE
  #is cycle on of the strange ones
  if(cycle %in% strange_cycle){ #2b0, 3b0,4b0 or 4u 
    
    if(length(grep("b", cycle))==1){
      zero_background=TRUE
      cycle=as.numeric(strsplit(cycle, "b")[[1]][1])
      
    }else{
      #4u
      cycle=as.numeric(strsplit(cycle, "u")[[1]][1])
    }
    
  }else{
    cycle=as.numeric(cycle)
  }
  
  
  df %>% filter(symbol==.env$symbol & experiment==.env$experiment)
  
  metadata_Yin2017$SELEX_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                            df %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% 
                                              select(filename)
                                            )
  #metadata_Yin2017$SELEX_ZeroCycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",df %>% filter(symbol=="ZeroCycle" & ligand==.env$ligand) %>% select(filename))
  if(zero_background!=TRUE){
  metadata_Yin2017$SELEX_previous_cycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                                           df %>% filter(symbol==.env$symbol & experiment==.env$experiment & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle-1) %>% 
                                                             select(filename)
                                                           )
  }else{
    
    metadata_Yin2017$SELEX_previous_cycle_filename[i]=paste0(data_path,"input_libraries", "/submitted_ftp/",
                                                             df_input %>% filter(ligand==.env$ligand) %>% 
                                                               select(filename)
    )
    
    
  }
}


lookup=c(CSC_SELEX_filename="SELEX_filename",            
         CSC_SELEX_previous_cycle_filename="SELEX_previous_cycle_filename",
         CSC_SELEX_ZeroCycle_filename="SELEX_ZeroCycle_filename")

metadata_Yin2017 <- metadata_Yin2017 %>%
  rename(all_of(lookup) )



metadata_Yin2017$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp/", 
                                             "https://a3s.fi/Yin2017/", 
                                             metadata_Yin2017$CSC_SELEX_filename)


metadata_Yin2017$Allas_SELEX_previous_cycle_filename=gsub("/scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp/", 
                                             "https://a3s.fi/Yin2017/", 
                                             metadata_Yin2017$CSC_SELEX_previous_cycle_filename)
metadata_Yin2017$Allas_SELEX_ZeroCycle_filename=gsub("/scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp/", 
                                                      "https://a3s.fi/Yin2017/", 
                                                      metadata_Yin2017$CSC_SELEX_ZeroCycle_filename)


missing_ind=which(str_detect(metadata_Yin2017$CSC_SELEX_filename, "/character"))
temp=metadata_Yin2017[missing_ind,]
metadata_Yin2017$CSC_SELEX_filename[missing_ind]=NA
metadata_Yin2017$Allas_SELEX_filename[missing_ind]=NA

missing_ind=which(str_detect(metadata_Yin2017$CSC_SELEX_previous_cycle_filename, "/character"))
temp=metadata_Yin2017[missing_ind,]
metadata_Yin2017$CSC_SELEX_previous_cycle_filename[missing_ind]=NA
metadata_Yin2017$Allas_SELEX_previous_cycle_filename[missing_ind]=NA



#included_ind=which(!str_detect(metadata_Yin2017$CSC_SELEX_previous_cycle_filename, "/character"))
#metadata_Yin2017$cycle[included_ind] %>% table()



#missing_ind=which(str_detect(metadata_Yin2017$CSC_SELEX_ZeroCycle_filename, "/character"))
#metadata_Yin2017$CSC_SELEX_ZeroCycle_filename[missing_ind]=NA
#metadata_Yin2017$Allas_SELEX_ZeroCycle_filename[missing_ind]=NA

metadata_Yin2017 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)

#"ZSCAN5A_HT-SELEX_TCGCCC40NCAT_KR_NYGTCCCYCCCCAAANMN_2_2" "ZSCAN31_Methyl-HT-SELEX_TAGTCT40NGCA_KV_GCATAACYGCCCYGCKKCN_2_4"
#"ZSCAN31_HT-SELEX_TGGAGA40NCCA_KV_GCATAACKGCCCTGCKKCN_2_4" 

tmp = metadata_Yin2017 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle, seed,Allas_SELEX_filename, 
         Allas_SELEX_previous_cycle_filename)

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Yin2017_Allas.tsv",
            delim="\t")


tmp = metadata_Yin2017 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,seed,CSC_SELEX_filename, 
         CSC_SELEX_previous_cycle_filename)


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Yin2017_CSC.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Yin2017$CSC_SELEX_filename %>% unique() %>% length() #962

metadata_Yin2017 %>% nrow() #1161


exp_multip=names(which(metadata_Yin2017$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Yin2017 %>% filter(CSC_SELEX_filename %in% exp_multip) #358



tmp$CSC_SELEX_filename %>% table() %>% table()
#2   3   4 
#133  20   8 

which( (tmp$CSC_SELEX_filename %>% table() )==4)
