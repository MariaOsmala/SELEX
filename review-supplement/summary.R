
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(httr2)
library(tibble)
library(httr)

setwd("~/projects/SELEX")
rm(list=ls())
file="/Users/osmalama/projects/TFBS/Data/SELEX-motif-collection/Supplementary_Table_1_Submission_July2026.tsv"

metadata=read_delim(file)

names(metadata)

summary(metadata)


metadata$ID %>%
  unique() %>%
  length() # 3933


metadata$kmer_scoring_exists %>% table(useNA="always")
#FALSE  TRUE  <NA> 
#  40  3595   298
is.na(metadata$N_sig) %>% table()
is.na(metadata$N_bg) %>% table()
is.na(metadata$lambda) %>% table()
#FALSE  TRUE 
#3596   337 

metadata %>% filter(is.na(lambda)) %>% select(experiment) %>% table()
#CAP-SELEX        HT-SELEX Methyl-HT-SELEX 
#2              38             297 




metadata %>%
  filter(experiment!="CAP-SELEX" & representative == "yes") %>%
  nrow() # 439

metadata %>%
  filter(experiment=="CAP-SELEX" & representative == "yes") %>%
  nrow() # 793



metadata$ID %>%
  unique() %>%
  length() # 3933

metadata$type %>%
  table(useNA = "always") %>%
  as.data.frame()
#. Freq
# 1     composite 1441
# 2       dimeric  822
# 3     monomeric 1067
# 4       spacing  457
# 5 ssDNA binding   74
# 6    tetrameric   27
# 7      trimeric   19
# 8       unknown   26
# 9          <NA>    0


metadata %>%
  filter(representative == "yes") %>%
  select(type) %>%
  table(useNA = "always") %>%
  as.data.frame()
# type Freq
# 1     composite  485
# 2       dimeric  186
# 3     monomeric  212
# 4       spacing  308
# 5 ssDNA binding   18
# 6    tetrameric    9
# 7      trimeric    3
# 8       unknown   11
# 9          <NA>    0
metadata %>%
  filter(is.na(type)) %>%
  pull(ID)




metadata %>% filter(representative=="yes") %>% select(sorted_TF1_copies) %>% table(useNA="always")
#1  1?   2  2?   3   4 
#811  37 326  33   3   9 

metadata %>% filter(representative=="yes") %>% select(sorted_TF2_copies) %>% table(useNA="always")


# For how many heterodimeric proteins the number of copies is unknown

metadata %>%
  filter(experiment == "CAP-SELEX" & representative == "yes") %>%
  filter(grepl("\\?", TF1_copies) | grepl("\\?", TF2_copies)) %>%
  nrow() # 254 (129 representative)

metadata %>%
  filter(experiment == "CAP-SELEX") %>%
  nrow() # 1896

metadata %>%
  filter(experiment == "CAP-SELEX" & representative == "yes") %>%
  nrow() # 793
# Are there some cases with more than 4 proteins, no

table(rowSums(data.frame(X1 = as.numeric(gsub("\\?", "", metadata$TF1_copies)), X2 = as.numeric(gsub("\\?", "", metadata$TF2_copies))), na.rm = TRUE))
# 0    1    2    3    4 
# 28 1141 2133  568   63 


apply(metadata[which(metadata$experiment == "CAP-SELEX"), c("TF1_copies", "TF2_copies")], 1, paste0, collapse = "_") %>% unique()





# Stats -------------------------------------------------------------------

# For how many we were able to find protein sequences

nrow(metadata) # 3933

(is.na(metadata %>% filter(representative == "yes") %>% pull(`protein sequence 1`))) %>% table()
#FALSE  TRUE 
#1062   170 

is.na(metadata$`protein sequence 1`) %>% table()
#FALSE  TRUE 
#3326   607 



metadata %>%
  pull(kmer_scoring_exists) %>%
  table()
#FALSE  TRUE 
#40  3595 



metadata %>%
  filter(!is.na(Human_Ensemble_ID)) %>%
  filter(Human_Ensemble_ID == `TF1 Ensembl_ID`) %>%
  nrow() # 1534

metadata %>%
  filter(!is.na(Human_Ensemble_ID)) %>%
  filter(Human_Ensemble_ID == `TF1_ensembl_gene_id`) %>%
  nrow() # 1509


# Remove Methyl-HT-SELEX motifs, ChIP-seq motifs, and ssDNA binding, and unknown motifs

data <- metadata %>% filter(experiment != "Methyl-HT-SELEX" & study != "Morgunova2015" ) # 3635

data <- data %>% mutate(
  SELEX_data_available = if_else(!is.na(CSC_SELEX_filename) & !is.na(CSC_SELEX_background_filename), "yes", "no")
)

data %>%
  filter(is.na(submitted_ftp_signal)) %>%
  nrow() # 30

data %>%
  filter(is.na(fastq_ftp_signal)) %>%
  nrow() # 39

data %>%
  filter(is.na(CSC_SELEX_filename)) %>%
  nrow() # 30

data %>%
  filter(is.na(submitted_ftp_background)) %>%
  nrow() # 26

data %>%
  filter(is.na(fastq_ftp_background)) %>%
  nrow() # 26

data %>%
  filter(is.na(CSC_SELEX_background_filename)) %>%
  nrow() # 26

data %>%
  pull(SELEX_data_available) %>%
  table()
#  no  yes
# 39 3596 

# kmer_scores_available

# protein_sequences_available
data <- data %>% mutate(
  protein_sequences_available = if_else(!is.na(TF1_protein_sequence) |
                                          (experiment == "CAP-SELEX" & (!is.na(TF1_protein_sequence) & !is.na(TF2_protein_sequence))), "yes", "no")
)
data %>%
  pull(protein_sequences_available) %>%
  table()
# no  yes
# 27 3608 

# Stoichiometry available

data$TF1_copies %>% table(useNA = "always")
data$TF2_copies %>% table(useNA = "always")


data <- data %>% mutate(
  stoichiometry_available = case_when(
    experiment != "CAP-SELEX" & grepl("^\\d+$", TF1_copies) ~ "yes",
    experiment == "CAP-SELEX" & grepl("^\\d+$", TF1_copies) & grepl("^\\d+$", TF2_copies) ~ "yes",
    TRUE ~ "no"
  )
)
data %>%
  pull(stoichiometry_available) %>%
  table()
#  no  yes
# 282 3353 
# stoichiometry_available,kmer_scores_available,protein_sequences_available,SELEX_data_available

# Draw Venn diagram of four sets: SELEX data available, k-mer scores available, protein sequences available, stoichiometry available

library(ggVennDiagram)
library(UpSetR)
library(ggplot2)

venn_sets <- list(
  "SELEX data\navailable" = data$ID[data$SELEX_data_available == "yes"],
  "k-mer scores\navailable" = data$ID[data$kmer_scoring_exists == TRUE],
  "Protein sequence(s)\navailable" = data$ID[data$protein_sequences_available == "yes"],
  "Stoichiometry\navailable" = data$ID[data$stoichiometry_available == "yes"]
)

# Venn diagram
venn_plot <- ggVennDiagram(venn_sets, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none")

ggsave(
  filename = "Figures/S2A_venn.pdf",
  plot = venn_plot,
  device = "pdf",
  width = 10,
  height = 8,
  units = "in"
)


# UpSet plot
upset_data <- data %>%
  transmute(
    `SELEX data available` = as.integer(SELEX_data_available == "yes"),
    `k-mer scores available` = as.integer(kmer_scoring_exists),
    `Protein sequence(s) available` = as.integer(protein_sequences_available == "yes"),
    `Stoichiometry available` = as.integer(stoichiometry_available == "yes")
  ) %>%
  as.data.frame()
pdf("Figures/S2A_upset.pdf", width = 10, height = 7, onefile = FALSE)
upset(
  upset_data,
  sets = c("SELEX data available", "k-mer scores available", "Protein sequence(s) available", "Stoichiometry available"),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  text.scale = 1.3
)
dev.off()

data <- data %>% filter(kmer_scoring_exists == TRUE & 
                          SELEX_data_available == "yes" & 
                          protein_sequences_available == "yes" & 
                          stoichiometry_available == "yes") # 3290

data <- data %>%
  mutate(
    CSC_SELEX_basename = basename(CSC_SELEX_filename),
    ENA_basename_signal = basename(submitted_ftp_signal)
  )

data %>%
  mutate(basename_signal = CSC_SELEX_basename == ENA_basename_signal) %>%
  select(basename_signal) %>%
  table() # all true

data <- data %>%
  mutate(
    CSC_SELEX_background_basename = basename(CSC_SELEX_background_filename),
    ENA_basename_background_signal = basename(submitted_ftp_background)
  )

data %>%
  mutate(basename_background = CSC_SELEX_background_basename == ENA_basename_background_signal) %>%
  select(basename_background) %>%
  table() # all true

#3290-13=3277

# There is some issue with this data
# 1 MGA_DLX2_CAP-SELEX_TGCGGT40NTCA_AAD_AGGTGNTAATTR_1_2u                  
# 2 POU2F1_ELK1_CAP-SELEX_TGACGA40NGCA_AS_NCCGGATATGCAN_1_2u             
# 3 POU2F1_ETV1_CAP-SELEX_TGCGAA40NAGC_AS_NCCGGATATGCAN_1_2u               
# 4 ERF_SREBF2_CAP-SELEX_TCTTTG40NACT_AAC_NNCACGTGACMGGAARNN_1_3u        
# 5 GCM1_ERG_CAP-SELEX_TAAGAA40NATA_AX_RTRYGGGCGGAARKN_1_3u                 
# 6 HOXA3_PAX5_CAP-SELEX_TCTGTC40NCAA_AY_YNATTAGTCACGCWTSRNTR_2_3u          
# 7 GCM2_DLX2_CAP-SELEX_TGTGCT40NCGG_AAB_RTRCGGGNNNNNTAATTR_1_3u             
# 8 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_ATRCGGGYNNNNNYWTTGTNN_1_3u         
# 9 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNNNNNYWTTGTNN_1_3u        
# 10 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNNRNACAAWN_1_3u           
# 11 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNRNACAAWN_1_3u             
# 12 HOXB7_HT-SELEX_TTATTT40NCGT_KV_NGTAATTANN_1_4u                           
# 13 TBX18_HT-SELEX_TATTTG40NCCA_KP_NRAGGTGTGAAN_1_4u                    



# What columns to keep in the table sent to the community?

cols_to_keep <- c(
  "motif_ID",
  "TF1", "TF1_ensembl_gene_id", "organism1", "TF1_copies", "TF1_protein_sequence",
  "TF2", "TF2_ensembl_gene_id", "organism2", "TF2_copies", "TF2_protein_sequence",
  "protein sequence 1", "protein sequence 2", "protein sequence 3", "protein sequence 4",
  "type",
  "Lambert2018_families",
  "study", "PMID", "experiment", "ligand",
  "seed",
  "cycle", "cycle_background", "unique_background",
  "representative",
  "filename",
  "run_accession_signal", "study_accession_signal",
  "secondary_study_accession_signal", "sample_accession_signal", "secondary_sample_accession_signal",
  "experiment_accession_signal", "submission_accession_signal", "read_count_signal",
  "base_count_signal", "submitted_ftp_signal", "fastq_ftp_signal", "sra_ftp_signal",
  "run_accession_background", "study_accession_background",
  "secondary_study_accession_background", "sample_accession_background", "secondary_sample_accession_background",
  "experiment_accession_background", "submission_accession_background", "read_count_background",
  "base_count_background", "submitted_ftp_background", "fastq_ftp_background", "sra_ftp_background",
  "TF1_dbd_matches_claimed_gene", "TF2_dbd_matches_claimed_gene",
  "kmer_scoring_exists", "SELEX_data_available",
  "protein_sequences_available", "stoichiometry_available"
)

data <- data %>% select(cols_to_keep)

# rename some column, use dplyr::rename

data <- data %>% rename(motif_filename = filename)

data <- data %>% rename(unique = unique_background)

data %>%
  select(TF1_copies) %>%
  table()
data %>%
  select(TF2_copies) %>%
  table()
data %>%
  select(type) %>%
  table()

data$TF1_dbd_matches_claimed_gene %>% table(useNA = "always")
data %>%
  filter(experiment == "CAP-SELEX") %>%
  select(TF2_dbd_matches_claimed_gene) %>%
  table(useNA = "always")
data %>%
  filter(experiment == "CAP-SELEX") %>%
  select(TF1_dbd_matches_claimed_gene) %>%
  table(useNA = "always")

# Remove columns

data <- data %>% select(-TF1_dbd_matches_claimed_gene, -TF2_dbd_matches_claimed_gene, -kmer_scoring_exists, -SELEX_data_available, -protein_sequences_available, -stoichiometry_available)

sum(data$read_count_signal, na.rm = TRUE) + sum(data$read_count_background, na.rm = TRUE) # 1 815 715 475 1.8 billion reads

data %>%
  filter(read_count_signal == 0) %>%
  pull(motif_ID)

# split Lambert2018_families into TF1_Lambert2018_families and TF2_Lambert2018_families

data <- data %>%
  separate(Lambert2018_families, into = c("TF1_family", "TF2_family"), sep = "_")

# TF1_family
# TF2_family

# Relocate TF1_family and TF2_family after TF1 and TF2

data <- data %>%
  relocate(TF1_family, .after = TF1) %>%
  relocate(TF2_family, .after = TF2)

data$motif_filename <- gsub("../../", "", data$motif_filename)

# write_delim(data,
#   "motif-metadata/S2A_curated_01042026.tsv",
#   delim = "\t"
# )

# write_delim(data,
#             "motif-metadata/S2A_curated_01072026.tsv",
#             delim = "\t"
# )

# 
# write_delim(metadata,
#   "motif-metadata/motifs_with_SELEX_signal_and_background_and_proteins_01042026.tsv",
#   delim = "\t"
# )

# write_delim(metadata,
#             "motif-metadata/motifs_with_SELEX_signal_and_background_and_proteins_01072026.tsv",
#             delim = "\t"
# )
