
library("readr")
library("tidyverse")
rm(list=ls())
df_motif_info=read_delim("~/projects/TFBS/Data/SELEX-motif-collection/metadata_final.tsv", "\t")



# Jolma2013 ---------------------------------------------------------------

df_Jolma2013 <- df_motif_info %>% filter(study=="Jolma2013")

# Replace instances of FL in column clone by "full"

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`clone` = ifelse(`clone` == "FL", "full", `clone`))

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`symbol` = ifelse(`clone` == "mouse DBD_mutant_DBD", "Egr1_E410D_FARSDERtoFARSDDR", `symbol`))

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`clone` = ifelse(`clone` == "mouse DBD_mutant_DBD", "DBD", `clone`))

df_Jolma2013=df_Jolma2013 %>% 
  mutate(symbol_clone = paste0(symbol,"_", clone)) %>% 
  relocate(symbol_clone, .after = clone)


## Read protein sequence info ----------------------------------------------

# Jolma 2013: Table S1 Protein and DNA sequence information, related to Figure 1

# A: Sequences of proteins used and their source; 
# the amino-acid sequence is based on end-sequencing of DNA and ENSEMBL. 
#
# Clones are either 
#   synthesized (Gene synthesis), 
#   PCR-cloned from the indicated source, 
#   or Gateway-cloned from the human ORFeome (hORFeome).

protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/1-s2.0-S0092867412014961-mmc1.xls", 
                                             skip = 3, 
                                             n_max=539)

#Correct some protein names
protein_sequences <- protein_sequences %>%
  mutate(`HNGC-name` = ifelse(`HNGC-name` == "Trp53", "Tp53", `HNGC-name`))

protein_sequences <- protein_sequences %>%
  mutate(`HNGC-name` = ifelse(`HNGC-name` == "Trp73", "Tp73", `HNGC-name`))

protein_sequences <- protein_sequences %>%
  mutate(symbol_clone = paste0(`HNGC-name`, "_", `clone type`) ) 

#Are the protein sequences unique for the same symbol_clone?
protein_sequences$symbol_clone %>% length() #538
protein_sequences$symbol_clone %>% unique() %>% length() #538
protein_sequences$`amino-acid sequence`%>% unique() %>% length() #534

df_Jolma2013_protein = df_Jolma2013 %>%
  left_join(protein_sequences, by="symbol_clone") 

#Did we find protein for all motifs YES!
df_Jolma2013_protein$`amino-acid sequence` %>% is.na() %>% table(useNA="always")


#Ensembl IDS no not always match and some missing
(df_Jolma2013_protein$Human_Ensemble_ID==df_Jolma2013_protein$`Ensembl id.`) %>% table(useNA="always")
#FALSE  TRUE  <NA> 
#  20   654   146 

#146 missing because it was not in Lambert table, are these mouse motifs. YES and some human motifs. 
df_Jolma2013_protein$Human_Ensemble_ID %>% is.na() %>% table()
#FALSE  TRUE  
#  674   146  
df_Jolma2013_protein %>% filter(is.na(Human_Ensemble_ID)) %>% select(organism) %>% table()
#Homo_sapiens Mus_musculus 
#13          133 

names(df_Jolma2013_protein)
#HNGC-name"            "Ensembl id."                
#"amino-acid sequence" "main structural class" "source organism" "clone type (DBD, full)"                 
# "production method"  "Clone source"               

df_Jolma2013_protein$`Classificaton in Vaquerizas et al., 2009`=NA

#Try to equalise the names 

df_Jolma2013_protein <- df_Jolma2013_protein %>%
  rename(`HNGC` = `HNGC-name`) 

df_Jolma2013_protein <- df_Jolma2013_protein %>%
  rename(`Ensembl_ID` = `Ensembl id.`)

df_Jolma2013_protein <- df_Jolma2013_protein %>%
  rename(`protein_sequence` = `amino-acid sequence`)

df_Jolma2013_protein <- df_Jolma2013_protein %>%
  rename(`DOMAIN` = `main structural class`)

Jolma2013_protein_sequences_unique=protein_sequences

#Try to add info to protein sequences (538 x 9)
Jolma2013_protein_sequences=protein_sequences %>% #820 x 25
  left_join(df_Jolma2013 %>% 
              select(ID,symbol_clone,study,experiment, family, Lambert2018_families, 
                     ligand, batch, seed, multinomial, cycle, representative, type,comment,
                     filename,Methyl.SELEX.Motif.Category, Human_Ensemble_ID
                
              ), by="symbol_clone") 

rm(protein_sequences)

Jolma2013_protein_sequences$`Classificaton in Vaquerizas et al., 2009`=NA

#Try to equalise the names 

Jolma2013_protein_sequences <- Jolma2013_protein_sequences %>%
  rename(`HNGC` = `HNGC-name`) 

Jolma2013_protein_sequences <- Jolma2013_protein_sequences %>%
  rename(`Ensembl_ID` = `Ensembl id.`)

Jolma2013_protein_sequences <- Jolma2013_protein_sequences %>%
  rename(`protein_sequence` = `amino-acid sequence`)

Jolma2013_protein_sequences <- Jolma2013_protein_sequences %>%
  rename(`DOMAIN` = `main structural class`)

Jolma2013_protein_sequences_unique <- Jolma2013_protein_sequences_unique %>%
  rename(`HNGC` = `HNGC-name`) 

Jolma2013_protein_sequences_unique <- Jolma2013_protein_sequences_unique %>%
  rename(`Ensembl_ID` = `Ensembl id.`)

Jolma2013_protein_sequences_unique <- Jolma2013_protein_sequences_unique %>%
  rename(`protein_sequence` = `amino-acid sequence`)

Jolma2013_protein_sequences_unique <- Jolma2013_protein_sequences_unique %>%
  rename(`DOMAIN` = `main structural class`)



#Jolma2013_protein_sequences %>% select(`DOMAIN`,family, Lambert2018_families) %>% View()

# Yin 2017 --------------------------------------------------------------

df_Yin2017 <- df_motif_info %>% filter(study=="Yin2017" & experiment=="HT-SELEX")

df_Yin2017=df_Yin2017 %>% 
  mutate(symbol_clone = paste0(symbol,"_", clone)) %>% 
  relocate(symbol_clone, .after = clone)

#Table S1 - Sequence information for Proteins (Section A, top) 	

#Section a		
#This table contains information about all of the TF constructs used in current study. 
#Columns from left to right contain:		
#Ensembl id.	Gene identifier	
#HNGC-name	HUGO Gene Nomenclature Committee accepted symbol for the TF	
#Construct type	Full length (FL) or extended DNA Binding Domain (eDBD)	
#main structural class	Main DNA binding domain in the construct	
#amino-acid sequence		
#Classification in Ref. Vaquerizas et al., 2009		


## Read protein sequence info ----------------------------------------------

protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/aaj2239_yin_sm_tables_s1-s6.xlsx", 
                                         sheet = "S1 sequence information", skip = 17, n_max= 1178-17)

protein_sequences <- protein_sequences %>%
  # Fix the column name typo
  rename(`eDBD sequence` = `eDBD seuqence`) %>% 
  separate_rows(`Clone used`, sep = ";\\s*")  # split into separate rows on ";"

protein_sequences <- protein_sequences %>%
  mutate(
    `eDBD sequence` = ifelse(`Clone used` == "eDBD", `eDBD sequence`, NA),
    `Full-length sequence` = ifelse(`Clone used` == "FL", `Full-length sequence`, NA)
  ) %>%
mutate(protein_sequence = coalesce(`eDBD sequence`, `Full-length sequence`)) %>%
  # Optionally remove the original two columns
  select(-`eDBD sequence`, -`Full-length sequence`)

protein_sequences <- protein_sequences %>%
  mutate(symbol_clone = paste0(`HNGC`, "_", `Clone used`) ) 

#Test for whitespaces in the symbol_clone column

protein_sequences$symbol_clone=str_replace_all(protein_sequences$symbol_clone, "\\s+", "")

#There is * is protein sequences? Remove it
grep("\\*", protein_sequences$protein_sequence) 

protein_sequences$protein_sequence=gsub("\\*","", protein_sequences$protein_sequence) 

#which(str_detect(df_Yin2017$symbol_clone, "\\s") )
#which(str_detect(protein_sequences$symbol_clone, "\\s") )



df_Yin2017_protein=df_Yin2017 %>%
  left_join(protein_sequences, by="symbol_clone", relationship = "many-to-many") 

test=df_Yin2017_protein[df_Yin2017_protein$protein_sequence %>% is.na() %>% which(),]

test %>% select(symbol_clone, Human_Ensemble_ID, Ensembl_ID)

#symbol_clone Human_Ensemble_ID Ensembl_ID     
#<chr>        <chr>             <chr>          
#1 SIX2_eDBD    ENSG00000170577   ENSG00000170577                           No eDBD protein, only full
#2 SPIB_FL      ENSG00000269404   ENSG00000142539    ENSG00000269404 SPIB   No full protein, only eDBD





#Maybe eDBD and full are the same?

#The following have different names in the protein_sequences table,
#the ensembl IDs are the same

#3 FOXO3_eDBD   ENSG00000118689   NA                 ENSG00000118689 FOXO3A
#4 FOXG1_eDBD   ENSG00000176165   NA                 ENSG00000176165 FOXG1B
#6 SKOR1_eDBD   ENSG00000188779   NA                 ENSG00000188779 LBXCOR1  
#7 ZSCAN5A_eDBD ENSG00000131848   NA                 ENSG00000131848 ZSCAN5    


missing_ind=df_Yin2017_protein$protein_sequence %>% is.na() %>% which()

df_Yin2017_protein$symbol[missing_ind]

lookup=c(FOXO3="FOXO3A", FOXG1="FOXG1B", SKOR1="LBXCOR1", ZSCAN5A="ZSCAN5")

for(i in 1:length(lookup)){
missing_ind=which(df_Yin2017_protein$symbol==names(lookup[i]))

#test=df_Yin2017_protein[missing_ind, ]

df_Yin2017_protein$symbol_clone[missing_ind]=gsub(names(lookup[i]),lookup[i],df_Yin2017_protein$symbol_clone[missing_ind])
df_Yin2017$symbol_clone[missing_ind]=gsub(names(lookup[i]),lookup[i],df_Yin2017$symbol_clone[missing_ind])

df_Yin2017_protein[missing_ind,]=df_Yin2017_protein[missing_ind,1:ncol(df_Yin2017)] %>%
  left_join(protein_sequences, by="symbol_clone", relationship = "many-to-many")

}

test=df_Yin2017_protein[df_Yin2017_protein$protein_sequence %>% is.na() %>% which(),]

test %>% select(symbol_clone, Human_Ensemble_ID, Ensembl_ID)

names(df_Yin2017_protein)
#"Ensembl_ID"  "HNGC"    "Clone used(eDBD   FL)"                              
# "DOMAIN"  "Classificaton in Vaquerizas et al., 2009" "protein_sequence"   

df_Yin2017_protein$`source organism`="human"

df_Yin2017_protein <- df_Yin2017_protein %>%
  rename(`clone type` = `Clone used`) 


df_Yin2017_protein$`production method`=NA
df_Yin2017_protein$`Clone source`=NA


#Is there some for which the protein info is missing: 

df_Yin2017_protein$protein_sequence %>% is.na() %>% table(useNA="always")
#FALSE  TRUE  <NA> 
#  865     2    0 

df_Yin2017_protein %>% filter(is.na(protein_sequence)) %>% select(ID, symbol_clone)
# SIX2_HT-SELEX_TTACTA40NAGG_KV_NCGTATCRYN_1_4    SIX2_eDBD   
# SPIB_HT-SELEX_TTTAGC40NGAC_KV_RAWWGMGGAAGTN_1_2 SPIB_FL   

#six_sequence="MSMLPTFGFTQEQVACVCEVLQQGGNIERLGRFLWSLPACEHLHKNESVLKAKAVVAFHRGNFRELYKILESHQFSPHNHAKLQQLWLKAHYIEAEKLRGRPLGAVGKYRVRRKFPLPRSIWDGEETSYCFKEKSRSVLREWYAHNPYPSPREKRELAEATGLTTTQVSNWFKNRRQRDRAAEAKERENNENSNSNSHNPLNGSGKSVLGSSEDEKTPSGTPDHSSSSPALLLSPPPPGLPSLHSLGHPPGPSAVPVPVPGGGGADPLQHHHGLQDSILNPMSANLVDLGS"
#nchar(six_sequence) #291

#protein_sequences %>% filter(HNGC %in% c("SIX2") & !is.na(protein_sequence)) %>% pull(protein_sequence) == six_sequence

#protein_sequences %>% filter(HNGC %in% c("SPIB") & !is.na(protein_sequence)) %>% pull(protein_sequence, symbol_clone)

#Try to add info to protein sequences (1534    7)

Yin2017_protein_sequences_unique=protein_sequences


Yin2017_protein_sequences=protein_sequences %>% #1728
  left_join(df_Yin2017 %>% 
              select(ID,symbol_clone,study,experiment, family, Lambert2018_families, 
                     ligand, batch, seed, multinomial, cycle, representative, type,comment,
                     filename,Methyl.SELEX.Motif.Category, Human_Ensemble_ID
                     
              ), by="symbol_clone", relationship = "many-to-many")





#How many motifs there are for each protein
#Yin2017_protein_sequences$number_of_motifs=sapply(Yin2017_protein_sequences$`df_Yin2017 %>% ...`, nrow)

Yin2017_protein_sequences %>% filter(!is.na(ID)) %>% nrow() #863

Yin2017_protein_sequences %>% filter(!is.na(ID)) %>% pull(symbol_clone) %>% unique() %>% length() #666

tmp=Yin2017_protein_sequences %>% filter(!is.na(ID)) %>% select(HNGC, `Clone used`) 

paste0(tmp$HNGC, "_", tmp$`Clone used`) %>% unique() %>% length() #666

# which(Yin2017_protein_sequences$symbol_clone %>% table(useNA="always") >1 ) %>% names()
#"GLIS2_eDBD" "MEF2B_eDBD" "PBX2_eDBD"  "SPIB_eDBD" 

#Are the protein sequences of these the same
#"GLIS2_eDBD(yes)" "MEF2B_eDBD (no, take the longer)" "PBX2_eDBD (yes)"  "SPIB_eDBD (yes)" 

#Yin2017_protein_sequences %>% filter(symbol_clone=="MEF2B_eDBD") %>% 
#  select(protein_sequence) %>% unique()

#Yin2017_protein_sequences$Ensembl_ID_alternative=NA

#index=which(Yin2017_protein_sequences$symbol_clone %in% c("GLIS2_eDBD","MEF2B_eDBD", "PBX2_eDBD",  "SPIB_eDBD"  ))

#keep_index=index[c(1,3,5,7)]
#remove_index=index[c(2,4,6,8)]

#Yin2017_protein_sequences$Ensembl_ID_alternative[keep_index]=Yin2017_protein_sequences$Ensembl_ID[remove_index]

#Yin2017_protein_sequences=Yin2017_protein_sequences[-remove_index,]

Yin2017_protein_sequences$`source organism`="human"

Yin2017_protein_sequences <- Yin2017_protein_sequences %>%
  rename(`clone type` = `Clone used`) 

Yin2017_protein_sequences$`production method`=NA
Yin2017_protein_sequences$`Clone source`=NA


Yin2017_protein_sequences_unique$`source organism`="human"

Yin2017_protein_sequences_unique <- Yin2017_protein_sequences_unique %>%
  rename(`clone type` = `Clone used`) 

Yin2017_protein_sequences_unique$`production method`=NA
Yin2017_protein_sequences_unique$`Clone source`=NA


# fam=sapply( Yin2017_protein_sequences$`df_Yin2017 %>% ...`, function(x) unique(x$Lambert2018_families))
# 
# fam[which(sapply(fam, length)==0)]=""
# Yin2017_protein_sequences$Lambert2018_families=unlist(fam)
# 
# fam=sapply( Yin2017_protein_sequences$`df_Yin2017 %>% ...`, function(x) unique(x$family))
# 
# fam[which(sapply(fam, length)==0)]=""
# Yin2017_protein_sequences$family=unlist(fam)


setdiff(names(Jolma2013_protein_sequences), names(Yin2017_protein_sequences))

Jolma2013_protein_sequences$symbol_clone %in% gsub("eDBD", "DBD", Yin2017_protein_sequences$symbol_clone) %>% table()
#FALSE  TRUE 
#375   445 

gsub("eDBD", "DBD", Yin2017_protein_sequences$symbol_clone)%in% Jolma2013_protein_sequences$symbol_clone  %>% table()
#FALSE  TRUE 
#1359   369 


#Are the protein sequences same in both studies

match_ind=match(Jolma2013_protein_sequences$symbol_clone, gsub("eDBD", "DBD", Yin2017_protein_sequences$symbol_clone))
na.ind=is.na(match_ind)

head(Jolma2013_protein_sequences$symbol_clone[-which(na.ind)],10)
head(gsub("eDBD", "DBD", Yin2017_protein_sequences$symbol_clone)[match_ind[-which(na.ind)]],10)

(Jolma2013_protein_sequences$protein_sequence[-which(na.ind)]==Yin2017_protein_sequences$protein_sequence[match_ind[-which(na.ind)]]) %>% table(useNA="always")
#FALSE  TRUE  <NA> 
#153   292   0

#Some are, some are not, maybe changes are not big ones

not_equal=which(!(Jolma2013_protein_sequences$protein_sequence[-which(na.ind)]==Yin2017_protein_sequences$protein_sequence[match_ind[-which(na.ind)]]))

Jolma2013_protein_sequences[ which(!na.ind)[not_equal[1]], c("symbol_clone", "protein_sequence")]
Yin2017_protein_sequences[match_ind[ which(!na.ind)][ not_equal[1]  ], c("symbol_clone", "protein_sequence")]

##Add Vaquerizas classification to Jolma2013 motifs

Vaquerizas_info <- Yin2017_protein_sequences %>%
  group_by(HNGC) %>%
  summarise(Va = unique(`Classificaton in Vaquerizas et al., 2009`), .groups = "drop")

match_ind=match(Jolma2013_protein_sequences$HNGC, Vaquerizas_info$HNGC) 
na.ind=is.na(match_ind)

head(Jolma2013_protein_sequences$HNGC[-which(na.ind)],50)
head(Vaquerizas_info$HNGC[match_ind[-which(na.ind)]],50)

Jolma2013_protein_sequences$`Classificaton in Vaquerizas et al., 2009`[-which(na.ind)] <- Vaquerizas_info$Va[match_ind[-which(na.ind)]]


# 
Jolma2013_protein_sequences$ID %>% unique() %>% length() #820

nrow(Yin2017_protein_sequences) #1728

Yin2017_protein_sequences$ID %>% unique() %>% length() #861
Yin2017_protein_sequences %>% filter(!is.na(ID)) %>% pull(ID) %>% unique() %>% length() #860, some missing here

#Add Vaqueriza classification to Jolma2013 motifs
df_Jolma2013_protein$protein_sequence %>% is.na() %>% table(useNA="always") #820 of 820

df_Yin2017_protein$protein_sequence %>% is.na() %>% table(useNA="always") #865. has protein, 2 missing

df_Jolma2013_protein %>% filter(symbol=="E2F2") %>% select(symbol_clone, protein_sequence)
df_Yin2017_protein %>% filter(symbol=="E2F2") %>% pull(symbol_clone,protein_sequence)

# Nitta2015 motifs --------------------------------------------------------

#Are the proteins in earlier studies
df_Nitta2015 <- df_motif_info %>% filter(study=="Nitta2015")

df_Nitta2015$symbol=str_replace_all(df_Nitta2015$symbol, "\\s+", "")

#Clone is missing

# Human HT-SELEX data were generated for PAX3, PGR, NR1I2, and NR1I3 using E. coli expressed DBDs,

df_Nitta2015$symbol
# "SREBF2"  "PAX4"    "PAX4"    "PAX3"    "PAX3"    "PGR"(not)     "ONECUT2" "ONECUT1" "NR1I2"   "NR1I3"

without_clone=c("SREBF2", "PAX4", "ONECUT2", "ONECUT1")

df_Nitta2015$symbol %in% Jolma2013_protein_sequences$HNGC
df_Nitta2015$symbol %in% Yin2017_protein_sequences$HNGC

Jolma2013_protein_sequences_unique %>% filter(HNGC %in% without_clone) %>% select(HNGC, protein_sequence, "clone type")
Yin2017_protein_sequences_unique %>% filter(HNGC %in% without_clone) %>% select(HNGC, protein_sequence, "clone type")

"ENSG00000082175" %in% Jolma2013_protein_sequences$Ensembl_ID
"ENSG00000082175" %in% Yin2017_protein_sequences$Ensembl_ID

# Jolma2015 ---------------------------------------------------------------

#Clone info missing?

#symbol(s)	HGNC symbol(s) (gene names); 
# Either for reference models that were generated for individual TFs in this study or 

#for the co-operative pairs, in which case the TF names are separated with an underscore "_". 
# First of the TFs is TF1 that was expressed as a SBP tagged clone and was used in first of the affinity separations 
# and the second one is the 3xFLAG tagged TF2	

df_Jolma2015_monomers <- df_motif_info %>% filter(study=="Jolma2015", experiment=="HT-SELEX") #31
df_Jolma2015_heterodimers <- df_motif_info %>% filter(study=="Jolma2015", experiment=="CAP-SELEX") #562

# Section a
# This table contains information about the TF1 and TF2 constructs that were used to run the CAP-SELEX.
# Columns from left to right contain:
# HNGC-name
# Ensembl id.
# amino-acid sequence
# main structural class
# 6-mer
# Clone source
# Construct type
# Protein amount
# HT-SELEX validation
# Activity in CAP-SELEX
# Included domain Ids
# Included domain description
# Excluded domain Ids
# Excluded domain description

protein_sequences=read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/41586_2015_BFnature15518_MOESM33_ESM.xlsx", 
                             sheet = "S1 sequence information", skip = 24,n_max =237-24)

protein_sequences$`Construct type` %>% table()
#Crystal Full length protein used in HT-SELEX with mixed proteins 
#2                                                       11 
#TF1(SBP)                                              TF2(3xFLAG) 
#105                                                       94 

#Remove Crustal and "Full length protein used in HT-SELEX with mixed proteins"

protein_sequences=protein_sequences %>% 
  filter(`Construct type` %in% c("TF1(SBP)", "TF2(3xFLAG)"))


#Are the protein sequences equal for replicate TFs (TF1(SBP)  and TF2(3xFLAG))?
multiple=names(which(protein_sequences$`HNGC-name` %>% table() >1))

#The sequences are the same are differ by 1 or 2 amino acids, consider the shorter

for(mtp in multiple){
  #mtp=multiple[1]
  seq=protein_sequences %>% filter(`HNGC-name`==mtp) %>% pull(`amino-acid sequence`)
  x=nchar(seq)
  #print(abs(outer(x, x, "-")))
  y=combn(x, 2, function(v) abs(diff(v)) )
  x_no_zeros <- y[y != 0]
  if(is.na(mean(x_no_zeros))){
    #no difference in length,  does the sequences differ, No they are all sam
    print(mtp)
    #print(all(seq ==seq[1]))
    #take the first
    tmp=protein_sequences %>% filter(`HNGC-name`==mtp)
    tmp=tmp[1,]
    protein_sequences=protein_sequences %>% filter(`HNGC-name`!=mtp)
    protein_sequences=rbind(protein_sequences, tmp)
    
    
    
    
    
  }else{
    #print(mtp)
    #print(mean(x_no_zeros)) #if the mean is 0, then all sequences are the same length
    #select some
    
    tmp=protein_sequences %>% filter(`HNGC-name`==mtp) 
    tmp=tmp[which.min(nchar(tmp$`amino-acid sequence`)),]
    
    protein_sequences=protein_sequences %>% filter(`HNGC-name`!=mtp)
    protein_sequences=rbind(protein_sequences, tmp)
    
    
  }
  #print(mean(x_no_zeros))
}



protein_sequences$symbol=protein_sequences$`HNGC-name`


df_Jolma2015_monomers_protein=df_Jolma2015_monomers %>%
  left_join(protein_sequences, by="symbol") 

df_Jolma2015_heterodimers= df_Jolma2015_heterodimers %>%
  rename("symbol_both"="symbol")

df_Jolma2015_heterodimers=df_Jolma2015_heterodimers %>%
          separate(symbol_both, into = c("symbol", "TF2"), sep = "_", remove=FALSE)

df_Jolma2015_heterodimers_protein=df_Jolma2015_heterodimers %>%
  left_join(protein_sequences, by="symbol") 

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename(
    `TF1 HNGC-name`=`HNGC-name`,
    `TF1 Ensembl id.`=`Ensembl id.`,
    `TF1 amino-acid sequence`=`amino-acid sequence`,
    `TF1 main structural class` =`main structural class` ,
    `TF1 6-mer`=`6-mer`,                      
    `TF1 source organism`=`source organism`,
    `TF1 Clone source`=`Clone source`,
    `TF1 Construct type`=`Construct type`, 
    `TF1 Protein amount`=`Protein amount`, 
    `TF1 HT-SELEX validation`=`HT-SELEX validation`,        
    `TF1 Activity in CAP-SELEX`=`Activity in CAP-SELEX`,
    `TF1 Included domain Ids`=`Included domain Ids`,
    `TF1 Included domain description`=`Included domain description`,
    `TF1 Excluded domain Ids`=`Excluded domain Ids`,
    `TF1 Excluded domain description`=`Excluded domain description`
  )

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename("TF1"="symbol")

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename("symbol"="TF2")

df_Jolma2015_heterodimers_protein=df_Jolma2015_heterodimers_protein %>%
  left_join(protein_sequences, by="symbol") 

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename(
    `TF2 HNGC-name`=`HNGC-name`,
    `TF2 Ensembl id.`=`Ensembl id.`,
    `TF2 amino-acid sequence`=`amino-acid sequence`,
    `TF2 main structural class` =`main structural class` ,
    `TF2 6-mer`=`6-mer`,                      
    `TF2 source organism`=`source organism`,
    `TF2 Clone source`=`Clone source`,
    `TF2 Construct type`=`Construct type`, 
    `TF2 Protein amount`=`Protein amount`, 
    `TF2 HT-SELEX validation`=`HT-SELEX validation`,        
    `TF2 Activity in CAP-SELEX`=`Activity in CAP-SELEX`,
    `TF2 Included domain Ids`=`Included domain Ids`,
    `TF2 Included domain description`=`Included domain description`,
    `TF2 Excluded domain Ids`=`Excluded domain Ids`,
    `TF2 Excluded domain description`=`Excluded domain description`
  )

#Everything ok

Jolma2015_protein_sequences=protein_sequences

# Xie 2025 ----------------------------------------------------------------

df_Xie2025_monomers <- df_motif_info %>% filter(study=="Xie2025" & experiment=="HT-SELEX") #12

#Is HHEX the same as THHEX?
df_Xie2025_monomers[which(df_Xie2025_monomers$symbol == "THHEX"), "symbol"] <- "HHEX"


df_Xie2025_heterodimers <- df_motif_info %>% filter(study=="Xie2025" & experiment=="CAP-SELEX") #1336

#These are unique proteins
protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/Xie et al. 2025 Supplementary_Tables.xlsx", 
                                sheet = "Table S1TF and ligand sequences", 
                                skip = 1270)

protein_sequences$`Clone used`=gsub("eDBD;", "eDBD", protein_sequences$`Clone used`)

protein_sequences$`Clone used` %>% table()
#eDBD   FL 
#430    3 

protein_sequences$symbol=protein_sequences$HNGC

df_Xie2025_monomers_protein=df_Xie2025_monomers %>%
  left_join(protein_sequences, by="symbol", relationship = "many-to-many") 

df_Xie2025_heterodimers= df_Xie2025_heterodimers %>%
  rename("symbol_both"="symbol")

df_Xie2025_heterodimers=df_Xie2025_heterodimers %>%
  separate(symbol_both, into = c("symbol", "TF2"), sep = "_", remove=FALSE)

df_Xie2025_heterodimers_protein=df_Xie2025_heterodimers %>%
  left_join(protein_sequences, by="symbol") 

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename(
    `TF1 Ensembl_ID.`=`Ensembl_ID`,
    `TF1 HNGC`=HNGC,
    `TF1 Clone used`=`Clone used`,
    `TF1 DOMAIN`=DOMAIN,
    `TF1 eDBD sequence`=`eDBD seuqence`
    
  )

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename("TF1"="symbol")

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename("symbol"="TF2")

df_Xie2025_heterodimers_protein=df_Xie2025_heterodimers_protein %>%
  left_join(protein_sequences, by="symbol") 

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename(
    `TF2 Ensembl_ID.`=`Ensembl_ID`,
    `TF2 HNGC`=HNGC,
    `TF2 Clone used`=`Clone used`,
    `TF2 DOMAIN`=DOMAIN,
    `TF2 eDBD sequence`=`eDBD seuqence`
    
  )

df_Xie2025_heterodimers_protein$`TF2 HNGC` %>% table(useNA="always")  #218
df_Xie2025_heterodimers_protein$`TF1 HNGC` %>% table(useNA="always")  #218

Xie2025_protein_sequences=protein_sequences


save(df_Jolma2013_protein,
       Jolma2013_protein_sequences,
       Jolma2013_protein_sequences_unique,
       df_Yin2017_protein,
       Yin2017_protein_sequences,
       Yin2017_protein_sequences_unique,
       df_Jolma2015_monomers_protein,
       df_Jolma2015_heterodimers_protein,
       Jolma2015_protein_sequences,
       df_Xie2025_monomers_protein,
       df_Xie2025_heterodimers_protein,
       Xie2025_protein_sequences, 
     file="~/projects/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData")


