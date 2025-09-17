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

#Nitta2015: AR_TCGCAA40NTGT_1.fastq.gz
#Morgunova2015: TFDP1.fastq.gz


metadata_Nitta2015 = metadata %>% filter(study=="Nitta2015")
study="Nitta2015"

filereport=read_delim("/scratch/project_2013895/SELEX/data/Nitta2015/filereport_read_run_PRJEB7373_tsv.txt")

filereport=filereport %>% filter(scientific_name=="Homo sapiens")



metadata_Nitta2015$SELEX_filename=NA
metadata_Nitta2015$SELEX_background_filename=NA

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
  
  metadata_Nitta2015$SELEX_filename[i]=paste0(data_path,study, "/submitted_ftp/", 
                                              filereport %>% filter(library_name==.env$symbol & sample_alias==.env$sample_alias ) %>% select(sample_alias), 
                                              ".fastq.gz")
  metadata_Nitta2015$SELEX_background_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                                         filereport %>% filter(library_name==.env$symbol & sample_alias==.env$sample_alias_background ) %>% select(sample_alias), 
                                                         ".fastq.gz")
}





lookup=c(CSC_SELEX_filename="SELEX_filename",    
         CSC_SELEX_background_filename="SELEX_background_filename")


metadata_Nitta2015 <- metadata_Nitta2015 %>%
  rename(all_of(lookup) )

metadata_Nitta2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp/", 
                                             "https://a3s.fi/Nitta2015/", 
                                             metadata_Nitta2015$CSC_SELEX_filename)


metadata_Nitta2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp/", 
                                             "https://a3s.fi/Nitta2015/", 
                                             metadata_Nitta2015$CSC_SELEX_background_filename)


missing_ind=which(str_detect(metadata_Nitta2015$CSC_SELEX_filename, "/character"))

temp=metadata_Nitta2015[missing_ind,]
temp$ID

# SREBF2\_HT-SELEX\_TGAGAT20NGA\_Y\_NRTCACGCCAYN\_1\_3
# ONECUT2\_HT-SELEX\_TGGGCG30NCGT\_AH\_NNCGATCRATAWNN\_1\_2
# ONECUT1\_HT-SELEX\_TAGCTC20NTCT\_Y\_CRATCRATAWN\_1\_3  

metadata_Nitta2015$CSC_SELEX_filename[missing_ind]=NA
metadata_Nitta2015$Allas_SELEX_filename[missing_ind]=NA


missing_ind=which(str_detect(metadata_Nitta2015$CSC_SELEX_background_filename, "/character"))
metadata_Nitta2015$CSC_SELEX_background_filename[missing_ind]=NA
metadata_Nitta2015$Allas_SELEX_background_filename[missing_ind]=NA

#background also missing

metadata_Nitta2015 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)
metadata_Nitta2015 %>% filter(is.na(CSC_SELEX_background_filename)) %>% pull(ID)


tmp = metadata_Nitta2015 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,cycle_background, unique_background,seed,CSC_SELEX_filename, 
         CSC_SELEX_background_filename, Allas_SELEX_filename, 
         Allas_SELEX_background_filename)


write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Nitta2015.tsv",
            delim="\t")


#Are the SELEX data listed several times

metadata_Nitta2015$CSC_SELEX_filename %>% unique() %>% length() #567

metadata_Nitta2015 %>% nrow() #820


exp_multip=names(which(metadata_Nitta2015$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Nitta2015 %>% filter(CSC_SELEX_filename %in% exp_multip) #175


tmp$CSC_SELEX_filename %>% table() %>% table()

