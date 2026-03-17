library("readr")
library("tidyverse")

rm(list=ls())

# Read metadata for HT-SELEX and CAP-SELEX motifs -------------------------

#motifs_with_SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background.tsv")
motifs_with_SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_06102025.tsv") #3594

#Are there duplicates

motifs_with_SELEX_data$ID %>% unique() %>% length() #3594


# For which motifs the score computations were successfull ----------------
lambda_info=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_lambda_info.tsv") #3594
lambda_info$ID %>% unique() %>% length() #3635

motifs_with_SELEX_data=motifs_with_SELEX_data %>% left_join(lambda_info %>% select(ID, lambda, N_sig, N_bg), by="ID")

is.na(motifs_with_SELEX_data$lambda) %>% table()
#FALSE  TRUE 
#3388   206 

is.na(motifs_with_SELEX_data$N_sig) %>% table() #124
is.na(motifs_with_SELEX_data$N_bg) %>% table() #124

# Read the type annotations -----------------------------------------------

# Monomer network

monomer_nodes=read_delim("/projappl/project_2013895/SELEX/cytoscape/monomer-network/node_table_24022026.csv", delim=",")

monomer_nodes=monomer_nodes %>% filter(node_type=="motif")

monomer_nodes %>% filter(representative=="YES") %>% nrow() #439

# Heterodimer network

heterodimer_nodes=read_delim("/projappl/project_2013895/SELEX/cytoscape/heterodimer-network/node_table_24022026.csv", delim=",")

heterodimer_nodes=heterodimer_nodes %>% filter(node_type=="motif")

heterodimer_nodes %>% filter(representative=="YES") %>% nrow() #793

motifs_with_SELEX_data =motifs_with_SELEX_data %>% left_join(monomer_nodes %>% select(name, representative,type, type_review), by=join_by("ID"=="name"))

motifs_with_SELEX_data =motifs_with_SELEX_data %>% left_join(heterodimer_nodes %>% select(name, representative,type, type_review), by=join_by("ID"=="name"))

motifs_with_SELEX_data$type=motifs_with_SELEX_data$type.x
motifs_with_SELEX_data$type[is.na(motifs_with_SELEX_data$type)]=motifs_with_SELEX_data$type.y[is.na(motifs_with_SELEX_data$type)]
motifs_with_SELEX_data$type.x=NULL
motifs_with_SELEX_data$type.y=NULL

motifs_with_SELEX_data$representative=motifs_with_SELEX_data$representative.x
motifs_with_SELEX_data$representative[is.na(motifs_with_SELEX_data$representative)]=motifs_with_SELEX_data$representative.y[is.na(motifs_with_SELEX_data$representative)]
motifs_with_SELEX_data$representative.x=NULL
motifs_with_SELEX_data$representative.y=NULL

motifs_with_SELEX_data$type_review=motifs_with_SELEX_data$type_review.x
motifs_with_SELEX_data$type_review[is.na(motifs_with_SELEX_data$type_review)]=motifs_with_SELEX_data$type_review.y[is.na(motifs_with_SELEX_data$type_review)]
motifs_with_SELEX_data$type_review.x=NULL
motifs_with_SELEX_data$type_review.y=NULL

motifs_with_SELEX_data$ID %>% unique() %>% length() #3594

motifs_with_SELEX_data$type %>% table(useNA="always") %>% as.data.frame()
motifs_with_SELEX_data %>% filter(representative=="YES") %>% select(type) %>% table(useNA="always") %>% as.data.frame()
motifs_with_SELEX_data %>% filter(is.na(type)) %>% pull(ID)
motifs_with_SELEX_data$type[which(is.na(motifs_with_SELEX_data$type))]=motifs_with_SELEX_data$type_review[which(is.na(motifs_with_SELEX_data$type))]

motifs_with_SELEX_data$type %>% table(useNA="always") %>% as.data.frame()
# . Freq
# 1       composite 1437
# 2      composite?    1
# 3         dimeric  709
# 4        dimeric?   19
# 5       monomeric  832
# 6      monomeric?    7
# 7         spacing  450
# 8        spacing?    8
# 9   ssDNA binding   72
# 10 ssDNA binding?    2
# 11     tetrameric   13
# 12    tetrameric?   10
# 13       trimeric   12
# 14        unknown   22
# 15           <NA>    0

motifs_with_SELEX_data  %>% pull(type_review) %>% table(useNA="always") %>% as.data.frame()

# 1      composite 1443
# 2        dimeric  726
# 3      monomeric  833
# 4     monomeric?    2
# 5        spacing  453
# 6  ssDNA binding   72
# 7     tetrameric   17
# 8       trimeric   18
# 9        unknown   22
# 10          <NA>    8


motifs_with_SELEX_data %>% filter(is.na(type_review)) %>% select(ID, representative)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  separate(symbol,
           into = c("TF1", "TF2"),
           sep = "_",
           fill = "right",
           remove = FALSE) %>%
  mutate(
    TF1_sorted = if_else(is.na(TF2), TF1, pmin(TF1, TF2)),
    TF2_sorted = if_else(is.na(TF2), NA_character_, pmax(TF1, TF2))
  ) %>%
  relocate(TF1_sorted, TF2_sorted, .after = TF2)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  mutate(sorted_TF1_copies = case_when(
    type_review == "monomeric" ~ 1,
    type_review == "dimeric" ~ 2,
    type_review == "trimeric" ~ 3,
    type_review == "tetrameric" ~ 4,
    TRUE ~ NA_real_
  ))


motifs_with_SELEX_data$sorted_TF2_copies=NA
#Check that the TF1_sorted and TF2_sorted match between tables 

heterodimer_ids=motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX") %>% pull(ID)
row_numbers <- which(motifs_with_SELEX_data$ID %in% heterodimer_ids)

match_inds=match(heterodimer_ids, heterodimer_nodes$name)

head(heterodimer_ids)
head(heterodimer_nodes$name[match_inds])

table((motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX") %>% pull(TF1_sorted))==heterodimer_nodes$TF1_sorted[match_inds]) #All TRUE
table((motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX") %>% pull(TF2_sorted))==heterodimer_nodes$TF2_sorted[match_inds]) #All TRUE

motifs_with_SELEX_data$sorted_TF1_copies[row_numbers]=heterodimer_nodes$sorted_TF1_copies[match_inds]
motifs_with_SELEX_data$sorted_TF2_copies[row_numbers]=heterodimer_nodes$sorted_TF2_copies[match_inds]

motifs_with_SELEX_data=motifs_with_SELEX_data %>% 
  relocate(sorted_TF1_copies, sorted_TF2_copies, type_review,.after=TF2_sorted)


# Add TF1_copies and TF2_copies 

motifs_with_SELEX_data$TF1_copies=motifs_with_SELEX_data$sorted_TF1_copies
motifs_with_SELEX_data$TF2_copies=motifs_with_SELEX_data$sorted_TF2_copies

swap_ids=row_numbers[which(motifs_with_SELEX_data$TF1[row_numbers]!=motifs_with_SELEX_data$TF1_sorted[row_numbers])]

motifs_with_SELEX_data$TF1_copies[swap_ids]=motifs_with_SELEX_data$sorted_TF2_copies[swap_ids]
motifs_with_SELEX_data$TF2_copies[swap_ids]=motifs_with_SELEX_data$sorted_TF1_copies[swap_ids]

motifs_with_SELEX_data=motifs_with_SELEX_data %>% 
  relocate(TF1_copies, TF2_copies, type_review,.after=TF2)

#For how many heterodimeric proteins the number of copies is unknown

motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX" & representative=="YES") %>% filter(grepl("\\?", TF1_copies) | grepl("\\?", TF2_copies)) %>% nrow() #253 (129 representative)
motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX") %>% nrow() #1896
motifs_with_SELEX_data %>% filter(experiment=="CAP-SELEX" & representative=="YES") %>% nrow() #793
#Are there some cases with more than 4 proteins

table( rowSums( data.frame(X1=as.numeric(gsub("\\?", "", motifs_with_SELEX_data$TF1_copies)), X2=as.numeric(gsub("\\?", "", motifs_with_SELEX_data$TF2_copies)) ) , na.rm=TRUE) )
#0    1    2    3    4 
#106  833 2035  567   53 

# Load protein sequence data ----------------------------------------------

load("/scratch/project_2013895/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData")

motifs_with_SELEX_data$`protein sequence 1`=NA
motifs_with_SELEX_data$`protein sequence 2`=NA
motifs_with_SELEX_data$`protein sequence 3`=NA
motifs_with_SELEX_data$`protein sequence 4`=NA


# Jolma2013 ---------------------------------------------------------------

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(
    df_Jolma2013_protein %>%
      select(ID, protein_sequence),
    by = "ID"
  ) %>%
  mutate(
    n = case_when(
      type_review == "monomeric"  ~ 1,
      type_review == "dimeric"    ~ 2,
      type_review == "trimeric"   ~ 3,
      type_review == "tetrameric" ~ 4,
      TRUE ~ 0
    ),
    
    `protein sequence 1` = if_else(n >= 1 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 1`),
    `protein sequence 2` = if_else(n >= 2 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 2`),
    `protein sequence 3` = if_else(n >= 3 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 3`),
    `protein sequence 4` = if_else(n >= 4 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 4`)
  ) %>%
  select(-protein_sequence, -n)





# Yin2017 -----------------------------------------------------------------

#df_Yin2017_protein$ID %>% unique() %>% length() #867 -> 864

#test=df_Yin2017_protein %>% filter(ID %in% names(which((df_Yin2017_protein$ID %>% table())>1)))

remove_ind=which(df_Yin2017_protein$ID %in% names(which((df_Yin2017_protein$ID %>% table())>1)))

df_Yin2017_protein=df_Yin2017_protein[-remove_ind[c(2,4,6)],]

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(
    df_Yin2017_protein %>%
      select(ID, protein_sequence),
    by = "ID"
  ) %>%
  mutate(
    n = case_when(
      type_review == "monomeric"  ~ 1,
      type_review == "dimeric"    ~ 2,
      type_review == "trimeric"   ~ 3,
      type_review == "tetrameric" ~ 4,
      TRUE ~ 0
    ),
    
    `protein sequence 1` = if_else(n >= 1 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 1`),
    `protein sequence 2` = if_else(n >= 2 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 2`),
    `protein sequence 3` = if_else(n >= 3 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 3`),
    `protein sequence 4` = if_else(n >= 4 & !is.na(protein_sequence),
                                   protein_sequence, `protein sequence 4`)
  ) %>%
  select(-protein_sequence, -n)


# Jolma2015 ---------------------------------------------------------------


## Monomers ----------------------------------------------------------------


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(
    df_Jolma2015_monomers_protein %>%
      select(ID, `amino-acid sequence`),
    by = "ID"
  ) %>%
  mutate(
    n = recode(type_review,
               monomeric  = 1L,
               dimeric    = 2L,
               trimeric   = 3L,
               tetrameric = 4L,
               .default   = 0L),
    
    `protein sequence 1` = if_else(n >= 1 & !is.na(`amino-acid sequence`),
                                   `amino-acid sequence`, `protein sequence 1`),
    
    `protein sequence 2` = if_else(n >= 2 & !is.na(`amino-acid sequence`),
                                   `amino-acid sequence`, `protein sequence 2`),
    
    `protein sequence 3` = if_else(n >= 3 & !is.na(`amino-acid sequence`),
                                   `amino-acid sequence`, `protein sequence 3`),
    
    `protein sequence 4` = if_else(n >= 4 & !is.na(`amino-acid sequence`),
                                   `amino-acid sequence`, `protein sequence 4`)
  ) %>%
  select(-`amino-acid sequence`, -n)


# Xie2025 -----------------------------------------------------------------


## Monomers ----------------------------------------------------------------

#Did not find SELEX data for monomers, do this does nothing

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(
    df_Xie2025_monomers_protein %>%
      select(ID, `eDBD seuqence`),
    by = "ID"
  ) %>%
  mutate(
    n = recode(type_review,
               monomeric  = 1L,
               dimeric    = 2L,
               trimeric   = 3L,
               tetrameric = 4L,
               .default   = 0L),
    
    `protein sequence 1` = if_else(n >= 1 & !is.na(`eDBD seuqence`),
                                   `eDBD seuqence`, `protein sequence 1`),
    
    `protein sequence 2` = if_else(n >= 2 & !is.na(`eDBD seuqence`),
                                   `eDBD seuqence`, `protein sequence 2`),
    
    `protein sequence 3` = if_else(n >= 3 & !is.na(`eDBD seuqence`),
                                   `eDBD seuqence`, `protein sequence 3`),
    
    `protein sequence 4` = if_else(n >= 4 & !is.na(`eDBD seuqence`),
                                   `eDBD seuqence`, `protein sequence 4`)
  ) %>%
  select(-`eDBD seuqence`, -n)



save.image("/scratch/project_2013895/SELEX/RData/tmp.RData")
load("/scratch/project_2013895/SELEX/RData/tmp.RData")
# Heterodimers ------------------------------------------------------------

## Jolma2015 Heterodimers ------------------------------------------------------------

#Rename column

df_Jolma2015_heterodimers_protein=df_Jolma2015_heterodimers_protein %>% rename(TF2=symbol)

#Check that the TF1 and TF2 match

match_ind=match(df_Jolma2015_heterodimers_protein$ID,motifs_with_SELEX_data$ID)
na_ind=which(is.na(match_ind))

(df_Jolma2015_heterodimers_protein$TF1[-na_ind]==motifs_with_SELEX_data$TF1[match_ind[-na_ind]]) %>% table() #all true
(df_Jolma2015_heterodimers_protein$TF2[-na_ind]==motifs_with_SELEX_data$TF2[match_ind[-na_ind]]) %>% table() #all true


#What kind of cases there can be: 1-1, 1-2, 2-1, 2-2, 1-3 

apply(motifs_with_SELEX_data[which(motifs_with_SELEX_data$experiment=="CAP-SELEX"),c("TF1_copies","TF2_copies")],1,paste0, collapse="_") %>% unique()

motifs_with_SELEX_data=motifs_with_SELEX_data %>% left_join(df_Jolma2015_heterodimers_protein %>% select(ID, `TF1 amino-acid sequence`, `TF2 amino-acid sequence`))

## Xie2025 Heterodimers ------------------------------------------------------------

View(df_Xie2025_heterodimers_protein)

df_Xie2025_heterodimers_protein=df_Xie2025_heterodimers_protein %>% rename(TF2=symbol)

#Check that the TF1 and TF2 match

match_ind=match(df_Xie2025_heterodimers_protein$ID,motifs_with_SELEX_data$ID)
na_ind=which(is.na(match_ind))

(df_Xie2025_heterodimers_protein$TF1[-na_ind]==motifs_with_SELEX_data$TF1[match_ind[-na_ind]]) %>% table() #all true
(df_Xie2025_heterodimers_protein$TF2[-na_ind]==motifs_with_SELEX_data$TF2[match_ind[-na_ind]]) %>% table() #all true


#What kind of cases there can be: 1-1, 1-2, 2-1, 2-2, 1-3 

apply(motifs_with_SELEX_data[which(motifs_with_SELEX_data$experiment=="CAP-SELEX"),c("TF1_copies","TF2_copies")],1,paste0, collapse="_") %>% unique()

motifs_with_SELEX_data=motifs_with_SELEX_data %>% left_join(df_Xie2025_heterodimers_protein %>% select(ID, `TF1 eDBD sequence`, `TF2 eDBD sequence`))

ind=which(!is.na(motifs_with_SELEX_data$`TF1 eDBD sequence`) | !is.na(motifs_with_SELEX_data$`TF2 eDBD sequence`))

motifs_with_SELEX_data$`TF1 amino-acid sequence`[ind]=motifs_with_SELEX_data$`TF1 eDBD sequence`[ind]
motifs_with_SELEX_data$`TF2 amino-acid sequence`[ind]=motifs_with_SELEX_data$`TF2 eDBD sequence`[ind]

motifs_with_SELEX_data$`TF1 eDBD sequence`=NULL
motifs_with_SELEX_data$`TF2 eDBD sequence`=NULL


motifs_with_SELEX_data=motifs_with_SELEX_data %>%
  rowwise() %>%
  mutate(
    valid = grepl("^\\d+$", TF1_copies) &
      grepl("^\\d+$", TF2_copies) &
      !is.na(`TF1 amino-acid sequence`) &
      !is.na(`TF2 amino-acid sequence`),
    
    proteins = list(
      if (valid) {
        x <- c(
          rep(`TF1 amino-acid sequence`, as.integer(TF1_copies)),
          rep(`TF2 amino-acid sequence`, as.integer(TF2_copies))
        )
        length(x) <- 4
        x
      } else {
        c(`protein sequence 1`,
          `protein sequence 2`,
          `protein sequence 3`,
          `protein sequence 4`)
      }
    ),
    
    `protein sequence 1` = proteins[1],
    `protein sequence 2` = proteins[2],
    `protein sequence 3` = proteins[3],
    `protein sequence 4` = proteins[4]
  ) %>%
  ungroup() %>%
  select(-valid, -proteins,
         -`TF1 amino-acid sequence`,
         -`TF2 amino-acid sequence`)



# Stats -------------------------------------------------------------------

#For how many we were able to find protein sequences

nrow(motifs_with_SELEX_data) #3594

(is.na(motifs_with_SELEX_data %>% filter(representative=="YES") %>% pull(`protein sequence 1`) )) %>% table()

is.na(motifs_with_SELEX_data$`protein sequence 1`) %>% table()
#FALSE  TRUE 
#3212   382

tmp=motifs_with_SELEX_data %>% select(ID, type_review, TF1, TF2, TF1_copies, TF2_copies, `protein sequence 1`, `protein sequence 2`, `protein sequence 3`, `protein sequence 4`)

#write_delim(motifs_with_SELEX_data, 
#            "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv", delim="\t")

# write_delim(motifs_with_SELEX_data, 
#             "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25022026.tsv", delim="\t")

write_delim(motifs_with_SELEX_data, 
            "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_04032026.tsv", delim="\t")



