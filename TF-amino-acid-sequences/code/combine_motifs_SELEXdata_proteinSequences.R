library("readr")
library("tidyverse")

rm(list=ls())
#only CAP-SELEX and HT-SELEX
motifs_with_SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background.tsv")

load("/scratch/project_2013895/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData")

motifs_with_SELEX_data$`protein sequence 1`=NA
motifs_with_SELEX_data$`protein sequence 2`=NA
motifs_with_SELEX_data$`protein sequence 3`=NA
motifs_with_SELEX_data$`protein sequence 4`=NA


# Jolma2013 ---------------------------------------------------------------

df_Jolma2013_protein$type %>% table(useNA="always")
#dimer of dimers               dimeric               monomer      monomer or dimer             monomeric 
#7                             434                   1            1                            355 
#monomeric or dimeric     putative multimer putatively multimeric              trimeric                  <NA> 
#1                        5                 9                                  7                     0 

df_Jolma2013_protein %>% filter(type %in% c("monomer", "monomeric", "monomer or dimer", "monomeric or dimeric")) 


motifs_with_SELEX_data=motifs_with_SELEX_data %>% #820 x 25
  left_join(df_Jolma2013_protein %>% 
              filter(type %in% c("monomer", "monomeric", "monomer or dimer", "monomeric or dimeric")) %>%
              select(ID, "protein_sequence", "type"),
            by="ID") %>%
mutate(`protein sequence 1` = protein_sequence) %>%
  select(-protein_sequence)

dimeric_data <- df_Jolma2013_protein %>%
  filter(type == "dimeric") %>%
  select(ID, protein_sequence)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(dimeric_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 2`), 
    type=if_else(!is.na(protein_sequence), "dimeric", `type`),
  ) %>%
  select(-protein_sequence)

trimeric_data <- df_Jolma2013_protein %>%
  filter(type == "trimeric") %>%
  select(ID, protein_sequence)


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(trimeric_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 2`), 
    `protein sequence 3` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 3`), 
    type=if_else(!is.na(protein_sequence), "trimeric", `type`),
  ) %>%
  select(-protein_sequence)

tetrameric_data <- df_Jolma2013_protein %>%
  filter(type == "dimer of dimers") %>%
  select(ID, protein_sequence)


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(tetrameric_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 2`), 
    `protein sequence 3` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 3`), 
    `protein sequence 4` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 4`), 
    type=if_else(!is.na(protein_sequence), "dimer of dimers", `type`),
  ) %>%
  select(-protein_sequence)

other_data <- df_Jolma2013_protein %>%
  filter(type %in% c("putative multimer", "putatively multimeric")) %>%
  select(ID, protein_sequence)

#Do not know what these are so add only once

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(other_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`)
  ) %>%
  select(-protein_sequence)

match_ind=match( motifs_with_SELEX_data$ID, other_data$ID)
notna_ind=which(!is.na(match_ind))

head(motifs_with_SELEX_data$ID[notna_ind],10)
head(other_data$ID[match_ind[notna_ind]],10)

motifs_with_SELEX_data$type[notna_ind]=(df_Jolma2013_protein %>%
  filter(type %in% c("putative multimer", "putatively multimeric")) %>% pull(type) )[match_ind[notna_ind]]


#P53 is tetramer, it is dimeric according to the table

tmp=motifs_with_SELEX_data %>% filter(study=="Jolma2013" & symbol=="Tp53")


# Yin2017 -----------------------------------------------------------------


Yin2017_data<- df_Yin2017_protein %>%
  select(ID, protein_sequence)

#Do not know what these are so add only once

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(Yin2017_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`)
  ) %>%
  select(-protein_sequence)


# Jolma2015 ---------------------------------------------------------------

#Did not find SELEX data for monomers
Jolma2015_data<- df_Jolma2015_monomers_protein %>%
  select(ID, `amino-acid sequence`)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(Jolma2015_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`amino-acid sequence`), `amino-acid sequence`, `protein sequence 1`)
  ) %>%
  select(-`amino-acid sequence`)


Jolma2015_data<- df_Jolma2015_heterodimers_protein %>%
  select(ID, `TF1 amino-acid sequence`, `TF2 amino-acid sequence`)


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(Jolma2015_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`TF1 amino-acid sequence`), `TF1 amino-acid sequence`, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(`TF2 amino-acid sequence`), `TF2 amino-acid sequence`, `protein sequence 2`)
  ) %>%
  select(-c(`TF1 amino-acid sequence`, `TF2 amino-acid sequence`))


# Xie2025 -----------------------------------------------------------------


#Did not find SELEX data for monomers
Xie2025_data<- df_Xie2025_monomers_protein %>%
  select(ID, `eDBD seuqence`)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(Xie2025_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`eDBD seuqence`), `eDBD seuqence`, `protein sequence 1`)
  ) %>%
  select(-`eDBD seuqence`)


Xie2025_data<- df_Xie2025_heterodimers_protein %>%
  select(ID, `TF1 eDBD sequence`, `TF2 eDBD sequence`)


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(Xie2025_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`TF1 eDBD sequence`), `TF1 eDBD sequence`, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(`TF2 eDBD sequence`), `TF2 eDBD sequence`, `protein sequence 2`)
  ) %>%
  select(-c(`TF1 eDBD sequence`, `TF2 eDBD sequence`))

#Add info for composite and spacing

match_ind=match( motifs_with_SELEX_data$ID, df_Xie2025_heterodimers_protein$ID)
notna_ind=which(!is.na(match_ind))

head(motifs_with_SELEX_data$ID[notna_ind],10)
head(df_Xie2025_heterodimers_protein$ID[match_ind[notna_ind]],10)

motifs_with_SELEX_data$type[notna_ind]=(df_Xie2025_heterodimers_protein %>% pull(type) )[match_ind[notna_ind]]


# Stats -------------------------------------------------------------------

#For how many we were able to find protein sequences

nrow(motifs_with_SELEX_data) #3576

is.na(motifs_with_SELEX_data$`protein sequence 1`) %>% table()

#3 missing which are there
#FALSE  TRUE 
#3566     10 

motifs_with_SELEX_data %>% filter( is.na(`protein sequence 1`)) %>% pull(ID)

write_delim(motifs_with_SELEX_data, 
            "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv", delim="\t")

