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

# Zero cycle background ---------------------------------------------------

input_library_filenames=dir(paste0(data_path,"input_libraries", "/submitted_ftp/"))

df_input <- as.data.frame(do.call(rbind, strsplit(input_library_filenames, "_")))

df_input$filename=input_library_filenames
df_input$V1=NULL
df_input$V3=NULL
df_input$V4=NULL
names(df_input)[1]="ligand"

#Remove those with no cycle reported
# cycle missing: BACH1_HT-SELEX_TTCCCC20NCCC_AL_ATGACTCAT_1_NA

metadata_Jolma2015=metadata_Jolma2015 %>% filter(!is.na(cycle))


metadata_Jolma2015$SELEX_filename=NA
metadata_Jolma2015$SELEX_background_cycle_filename=NA
metadata_Jolma2015$unique_background=FALSE
metadata_Jolma2015$cycle_background=NA

  
for(i in 1:nrow(metadata_Jolma2015)){
  #i=1
  print(i)
  symbol=metadata_Jolma2015$symbol[i]
  
  experiment=metadata_Jolma2015$experiment[i]
  ligand=metadata_Jolma2015$ligand[i]
  batch=metadata_Jolma2015$batch[i]
  cycle=metadata_Jolma2015$cycle[i]
  
  df %>% filter(symbol==.env$symbol)
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
  
  metadata_Jolma2015$SELEX_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                              df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle) %>% 
                                                select(filename)
                                              )
  if(cycle_background!=0){
  metadata_Jolma2015$SELEX_background_cycle_filename[i]=paste0(data_path,study, "/submitted_ftp/",
                                                             df %>% filter(symbol==.env$symbol & ligand==.env$ligand & batch==.env$batch & cycle==.env$cycle_background) %>% select(filename))
  }else{
    metadata_Jolma2015$SELEX_background_cycle_filename[i]=paste0(data_path,"input_libraries", "/submitted_ftp/",
                                                               df_input %>% filter(symbol==.env$symbol & ligand==.env$ligand) %>% select(filename))
  }
}







lookup=c(CSC_SELEX_filename="SELEX_filename",    
         CSC_SELEX_background_filename="SELEX_background_cycle_filename")


metadata_Jolma2015 <- metadata_Jolma2015 %>%
  rename(all_of(lookup) )

metadata_Jolma2015$Allas_SELEX_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2015/", 
                                             metadata_Jolma2015$CSC_SELEX_filename)


metadata_Jolma2015$Allas_SELEX_background_filename=gsub("/scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp/", 
                                             "https://a3s.fi/Jolma2015/", 
                                             metadata_Jolma2015$CSC_SELEX_background_filename)


missing_ind=which(str_detect(metadata_Jolma2015$CSC_SELEX_filename, "/character"))

tmp=metadata_Jolma2015[missing_ind,]

metadata_Jolma2015$CSC_SELEX_filename[missing_ind]=NA
metadata_Jolma2015$Allas_SELEX_filename[missing_ind]=NA

# CEBPG\_ATF4\_CAP-SELEX\_TGCGTC40NTTA\_AAB\_NNATGAYGCAAT\_1\_3b0 
# ELF2\_HT-SELEX\_TGCAAG20NAAC\_AL\_NAMCCGGAAGTR\_1\_2          
# ELF2\_HT-SELEX\_TGCAAG20NAAC\_AL\_NATGCGGAAGTR\_1\_2
# ETS2\_HT-SELEX\_TAAGTG40NGAA\_AR\_RCCGGAAGTG\_1\_2            
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNGCGGAAGTG\_1\_4            
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNGGAAGTGCTTCCNN\_2\_4      
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNYTTCCGGGAARNR\_1\_4       
# FOS\_HT-SELEX\_TGAACT40NAAG\_KR\_NGATGACGTCATCR\_2\_4         
# FOXA1\_HT-SELEX\_TTCTAA40NAAT\_KN\_TRNGTAAACA\_1\_3b1         
# GATA1\_HT-SELEX\_TGGGTA20NTGT\_AL\_AGATAAN\_1\_2              
# GLI3\_HT-SELEX\_TACCCG20NCCC\_AN\_NGACCACMCACGWNG\_2\_3      
# HES1\_HT-SELEX\_TCTTTC20NTTG\_AL\_GNCACGTGNC\_1\_3            
# HOXA3\_HT-SELEX\_TGTCGT40NGCG\_AL\_NSTAATTANN\_1\_3          
# HOXA4\_HT-SELEX\_TGACCT40NCCA\_AR\_RTMATTAN\_1\_4             
# HOXA6\_HT-SELEX\_TCGCCA20NGA\_AN\_SYMATTAN\_1\_3            
# HOXA7\_HT-SELEX\_TCGCGC20NGA\_AN\_NYMATTAN\_1\_3              
# HOXD4\_HT-SELEX\_TGGCCC40NCCT\_AR\_NNYMATTANN\_1\_4b0         
# IRF2\_HT-SELEX\_TTGCCC40NCTC\_AAF\_NAANCGAAASYR\_1\_3         
# JUN\_HT-SELEX\_TTAGCC20NTA\_AL\_ATGACGTCAT\_1\_3            
# MYCL2\_HT-SELEX\_TAGCCT40NCCT\_AR\_SCACGTGS\_1\_3             
# NR1D2\_HT-SELEX\_TGAATT40NTAA\_KR\_TRGGTYASTAGGTCA\_2\_3     
# NR1I2\_HT-SELEX\_TACTGG40NGGA\_KR\_RGTTCRNNNRGTTC\_1\_4       
# POU5F1\_HT-SELEX\_TTTAAG40NAGC\_AQ\_NATATGCTAATKN\_1\_3      
# POU5F1\_HT-SELEX\_TTTAAG40NAGC\_AQ\_WATGCGCATW\_1\_3          
# RORB\_HT-SELEX\_TTCGGG40NGAG\_KS\_AANTAGGTCAGTAGGTCA\_2\_4   
# RORB\_HT-SELEX\_TTCGGG40NGAG\_KS\_AWNTAGGTCATGACCTANWT\_2\_4  
# SOX17\_HT-SELEX\_TATGCT40NACT\_KO\_ACCGAACAAT\_1\_4b2        
# SOX6\_HT-SELEX\_TTCCAA20NACC\_AL\_CACCGAACAAT\_2\_3           
# TBX3\_HT-SELEX\_TAAGCC40NAGT\_AR\_AGGTGTNR\_1\_4            
# TCF15\_HT-SELEX\_TCTTAG40NATG\_AR\_NACAYATGNN\_1\_4           
# PAX3\_HT-SELEX\_TTAGGG20NGGA\_AL\_GTCACGCNNMATTAN\_1\_3 


missing_ind=which(str_detect(metadata_Jolma2015$CSC_SELEX_background_filename, "/character"))

tmp=metadata_Jolma2015[missing_ind,]

metadata_Jolma2015$CSC_SELEX_background_filename[missing_ind]=NA
metadata_Jolma2015$Allas_SELEX_background_filename[missing_ind]=NA

#background missing: 

# ELF2\_HT-SELEX\_TGCAAG20NAAC\_AL\_NAMCCGGAAGTR\_1\_2    
# ELF2\_HT-SELEX\_TGCAAG20NAAC\_AL\_NATGCGGAAGTR\_1\_2        
# ETS2\_HT-SELEX\_TAAGTG40NGAA\_AR\_RCCGGAAGTG\_1\_2          
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNGCGGAAGTG\_1\_4         
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNGGAAGTGCTTCCNN\_2\_4     
# ETV7\_HT-SELEX\_TCTGAT40NCTA\_AQ\_NNYTTCCGGGAARNR\_1\_4     
# FOS\_HT-SELEX\_TGAACT40NAAG\_KR\_NGATGACGTCATCR\_2\_4      
# FOXA1\_HT-SELEX\_TTCTAA40NAAT\_KN\_TRNGTAAACA\_1\_3b1       
# GATA1\_HT-SELEX\_TGGGTA20NTGT\_AL\_AGATAAN\_1\_2          
# GLI3\_HT-SELEX\_TACCCG20NCCC\_AN\_NGACCACMCACGWNG\_2\_3     
# HES1\_HT-SELEX\_TCTTTC20NTTG\_AL\_GNCACGTGNC\_1\_3        
# HOXA3\_HT-SELEX\_TGTCGT40NGCG\_AL\_NSTAATTANN\_1\_3         
# HOXA4\_HT-SELEX\_TGACCT40NCCA\_AR\_RTMATTAN\_1\_4          
# HOXA6\_HT-SELEX\_TCGCCA20NGA\_AN\_SYMATTAN\_1\_3            
# HOXA7\_HT-SELEX\_TCGCGC20NGA\_AN\_NYMATTAN\_1\_3           
# IRF2\_HT-SELEX\_TTGCCC40NCTC\_AAF\_NAANCGAAASYR\_1\_3       
# JUN\_HT-SELEX\_TTAGCC20NTA\_AL\_ATGACGTCAT\_1\_3           
# MYCL2\_HT-SELEX\_TAGCCT40NCCT\_AR\_SCACGTGS\_1\_3           
# NR1D2\_HT-SELEX\_TGAATT40NTAA\_KR\_TRGGTYASTAGGTCA\_2\_3   
# NR1I2\_HT-SELEX\_TACTGG40NGGA\_KR\_RGTTCRNNNRGTTC\_1\_4     
# POU5F1\_HT-SELEX\_TTTAAG40NAGC\_AQ\_NATATGCTAATKN\_1\_3     
# POU5F1\_HT-SELEX\_TTTAAG40NAGC\_AQ\_WATGCGCATW\_1\_3        
# RORB\_HT-SELEX\_TTCGGG40NGAG\_KS\_AANTAGGTCAGTAGGTCA\_2\_4 
# RORB\_HT-SELEX\_TTCGGG40NGAG\_KS\_AWNTAGGTCATGACCTANWT\_2\_4
# SOX17\_HT-SELEX\_TATGCT40NACT\_KO\_ACCGAACAAT\_1\_4b2      
# SOX6\_HT-SELEX\_TTCCAA20NACC\_AL\_CACCGAACAAT\_2\_3         
# TBX3\_HT-SELEX\_TAAGCC40NAGT\_AR\_AGGTGTNR\_1\_4           
# TCF15\_HT-SELEX\_TCTTAG40NATG\_AR\_NACAYATGNN\_1\_4         
# PAX3\_HT-SELEX\_TTAGGG20NGGA\_AL\_GTCACGCNNMATTAN\_1\_3     



#missing_ind=which(str_detect(metadata_Jolma2015$CSC_SELEX_ZeroCycle_filename, "/character"))
#metadata_Jolma2015$CSC_SELEX_ZeroCycle_filename[missing_ind]=NA
#metadata_Jolma2015$Allas_SELEX_ZeroCycle_filename[missing_ind]=NA

metadata_Jolma2015 %>% filter(is.na(CSC_SELEX_filename)) %>% pull(ID)


tmp = metadata_Jolma2015 %>% 
  select(ID,symbol,clone,Lambert2018_families,experiment,ligand, batch,cycle,cycle_background, unique_background, seed,CSC_SELEX_filename, 
         CSC_SELEX_background_filename, Allas_SELEX_filename, 
         Allas_SELEX_background_filename)

  
write_delim(tmp,"/projappl/project_2013895/motif_metadata/metadata_Jolma2015.tsv",
            delim="\t")



#Are the SELEX data listed several times

metadata_Jolma2015$CSC_SELEX_filename %>% unique() %>% length() #328

metadata_Jolma2015 %>% nrow() #592


exp_multip=names(which(metadata_Jolma2015$CSC_SELEX_filename %>% table() >1))

tmp=metadata_Jolma2015 %>% filter(CSC_SELEX_filename %in% exp_multip) #175

tmp$CSC_SELEX_filename %>% table() %>% table()
#2  3  4  5 
#87 37 15  7

which( (tmp$CSC_SELEX_filename %>% table() )==5)
