library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

rm(list=ls())


Lambert2018_supplement="Data/Lambert2018_supplementary_tables_S1_S4.xlsx"


Lambert2018 <- read_excel(Lambert2018_supplement,
                          sheet = "Table S1. Related to Figure 1B",
                          col_types="text",
                          skip = 1)

#Rename the fourth column

colnames(Lambert2018)[4]="Is TF?"

TFs=Lambert2018 %>% filter(`Is TF?`=="Yes") # 1639 TFs

#Could search for all TF protein info from supplementary tables
# Add proteins from Xie2025 paper

#Get the DBD and ensemblID info from Xie2025 supplement
#Section b
#This table contains information about the TF1 and TF2 constructs used to run the CAP-SELEX.
#Columns from left to right contain:
#Ensembl id.	Gene identifier
#HNGC-name	HUGO Gene Nomenclature Committee accepted symbol for the TF
#Construct type	TF1(His) or TF2(SBP) that are used in the first and second affinity purification steps, respectively
#main structural class	Main DNA binding domain in the construct
#amino-acid sequence	TF-specific part of the TF1 or TF2 construct
protein_info <- read_excel("~/projects/TFBS/Data/SELEX-motif-collection/41586_2025_8844_MOESM3_ESM.xlsx",
                           sheet="Supplementary Table S1", col_names=TRUE, skip=1269, n_max=434)

protein_info$HNGC %in% TFs$Name %>% table()

Xie2025_TFs=protein_info %>% filter(!(protein_info$HNGC %in% TFs$Name)) %>% dplyr::select(Ensembl_ID, HNGC, DOMAIN) #23
names(Xie2025_TFs)=c("ID", "Name", "DBD")
TFs=dplyr::bind_rows(TFs, Xie2025_TFs) # 1662

#Does TFs contain all proteins tested in SELEX

TFBS_path="/Users/osmalama/projects/TFBS/"
metadata <- read_tsv(paste0(TFBS_path,"Data/SELEX-motif-collection/metadata_final.tsv"), col_names=TRUE)

metadata_monomers=metadata %>% filter(experiment %in% c("HT-SELEX", "Methyl-HT-SELEX"))

# add new column after symbol
metadata_monomers$alter_symbol=""
metadata_monomers=metadata_monomers %>% relocate(alter_symbol, .after = symbol)

metadata_monomers$alter_symbol[which(metadata_monomers$symbol=="THHEX")]="THHEX"
metadata_monomers$symbol[which(metadata_monomers$symbol=="THHEX")]="HHEX"

metadata_heterodimers=metadata %>% filter(!(experiment %in% c("HT-SELEX", "Methyl-HT-SELEX")))
metadata_heterodimers= metadata_heterodimers %>% separate(symbol, into = c("TF1", "TF2"), sep = "_", remove=FALSE) %>%
  relocate(TF1, .after = symbol) %>% relocate( TF2, .after = TF1)

metadata_heterodimers= metadata_heterodimers %>%
  separate(Lambert2018_families,
    into = c("TF1_Lambert2018_families", "TF2_Lambert2018_families"), sep = "_", remove=FALSE) %>%
  relocate(TF1_Lambert2018_families, .after = Lambert2018_families) %>%
  relocate( TF2, .after = TF1_Lambert2018_families)

TFs_with_motifs=unique(c(metadata_monomers$symbol, metadata_heterodimers$TF1, metadata_heterodimers$TF2)) # 728

saveRDS(TFs_with_motifs, file="RData/TFs_with_motifs.RDS")

toupper(TFs_with_motifs) %in% TFs$Name %>% table()
#FALSE  TRUE
#13   715

missing_TFs=TFs_with_motifs[!(toupper(TFs_with_motifs) %in% TFs$Name)]

#What is the family of these missing TFs

missing_TFs_families=metadata_monomers %>% filter(symbol %in% unique(missing_TFs)) %>% dplyr::select(symbol, Lambert2018_families) %>% unique()

# Add new rows of TFs
names(missing_TFs_families)=c("Name", "DBD")

TFs=dplyr::bind_rows(TFs, missing_TFs_families) # 1674



#For how many TFs there is a motif in the SELEX data (monomers and/or monomultimers?)

#Homologous motifs?

(TFs %>% pull(Name) %in% toupper(metadata_monomers$symbol))%>% table()
#FALSE  TRUE
#1033   641

TFs$has_SELEX_motif=(TFs %>% pull(Name) %in% toupper(metadata_monomers$symbol))

TFs %>%
  count(has_SELEX_motif)



# Section a
# This table contains interaction information between TF1 (row, His-tagged) and TF2 (column, SBP-tagged).
# The numbers in the table describe the inteaction types:
# "0"	Tested TF pairs without cooperative signal
# "1"	The interaction between TFs results in composite motifs
# "2"	The interaction between TFs results in spacing and/or oreintaton preference
# N/A	TF pairs have not been tested in this study
#
# broad	A TF promiscuously interacts with TFs from the same or other structural family
# subfamily	A TF interacts with TFs from a distinct class of a TF family
# Paralog	A TF specifically interacts with paralogs from the same or other TF family
# Specific	A TF interacts with one or few proteins from  the same or other structural family

# 159 rows, 398 columns

# 158 TFs (His-tagged preys) tested against 419 TFs (SBP-tagged prey) (10% and 25% of human proteins, respectively
#In CAP-SELEX, 158 prey proteins were screened against a set of 376 bait proteins arranged on 384-well plates.

interaction_matrix <- read_excel("~/projects/TFBS/Data/SELEX-motif-collection/41586_2025_8844_MOESM4_ESM.xlsx",
                        sheet="Supplementary Table S2", col_names=TRUE, skip=31, n_max=159)


Preys=interaction_matrix$`Prey\\Bait` #rows #159
Baits=names(interaction_matrix)[-c(1,396,397,398,399)] #columns #394

interaction_matrix=as.data.frame(interaction_matrix[,-c(1,396,397,398,399)])
rownames(interaction_matrix)=Preys

#Are all Preys and Baits in TFs

unique(c(Preys[!(Preys %in% TFs$Name)], Baits[!(Baits %in% TFs$Name)])) #16

#"LRRN6D"
#"RFXDC2"

#which are missing
unique(c(Preys[!(Preys %in% TFs$Name)], Baits[!(Baits %in% TFs$Name)]))[!(unique(c(Preys[!(Preys %in% TFs$Name)], Baits[!(Baits %in% TFs$Name)])) %in% protein_info$HNGC)] #TRUE
#"LRRN6D"
#"RFXDC2"

extra=data.frame(ID=c("ENSG00000213171", "ENSG00000181827"), Name=c("LRRN6D", "RFXDC2"), DBD=c("LRR", "RFX"))
TFs=dplyr::bind_rows(TFs, extra)


