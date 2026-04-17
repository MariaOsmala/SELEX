library("readr")
library("tidyverse")

rm(list = ls())

# Read metadata for HT-SELEX and CAP-SELEX motifs -------------------------

# motifs_with_SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background.tsv")
# motifs_with_SELEX_data=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_06102025.tsv") #3594

motifs_with_SELEX_data <- read_delim("motif-metadata/metadata_final_with_SELEX_data_19032026.tsv") # 3933
# Are there duplicates

motifs_with_SELEX_data$ID %>%
  unique() %>%
  length() # 3933

# For which motifs the score computations were successfull ----------------
lambda_info <- read_delim("motif-metadata/metadata_final_with_SELEX_data_lambda_info.tsv") # 3594
lambda_info$ID %>%
  unique() %>%
  length() # 3635

motifs_with_SELEX_data <- motifs_with_SELEX_data %>% left_join(lambda_info %>% select(ID, lambda, N_sig, N_bg), by = "ID")

is.na(motifs_with_SELEX_data$lambda) %>% table()
# FALSE  TRUE
# 3388   545

is.na(motifs_with_SELEX_data$N_sig) %>% table() # 463
is.na(motifs_with_SELEX_data$N_bg) %>% table() # 463

# Read the type annotations -----------------------------------------------

# Monomer network

# monomer_nodes=read_delim("/projappl/project_2013895/SELEX/cytoscape/monomer-network/node_table_24022026.csv", delim=",")
monomer_nodes <- read_delim("../cytoscape/monomer_network/node_table_24022026.csv", delim = ",")

monomer_nodes <- monomer_nodes %>% filter(node_type == "motif")

monomer_nodes %>%
  filter(representative == "YES") %>%
  nrow() # 439

# Heterodimer network

# heterodimer_nodes=read_delim("/projappl/project_2013895/SELEX/cytoscape/heterodimer-network/node_table_24022026.csv", #delim=",")
heterodimer_nodes <- read_delim("../cytoscape/heterodimer_network/node_table_24022026.csv", delim = ",")

heterodimer_nodes <- heterodimer_nodes %>% filter(node_type == "motif")

heterodimer_nodes %>%
  filter(representative == "YES") %>%
  nrow() # 793

motifs_with_SELEX_data$type <- NULL

# Lambert2018_families are the same in both tables

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  left_join(
    monomer_nodes %>% select(name, type, type_review),
    by = join_by(ID == name)
  ) %>%
  left_join(
    heterodimer_nodes %>% select(name, type, type_review),
    by = join_by(ID == name),
    suffix = c("_mono", "_hetero")
  ) %>%
  mutate(
    type = coalesce(type_mono, type_hetero),
    type_review = coalesce(type_review_mono, type_review_hetero)
  ) %>%
  select(-type_mono, -type_hetero, -type_review_mono, -type_review_hetero)



motifs_with_SELEX_data$ID %>%
  unique() %>%
  length() # 3933

motifs_with_SELEX_data$type %>%
  table(useNA = "always") %>%
  as.data.frame()
motifs_with_SELEX_data %>%
  filter(representative == "YES") %>%
  select(type) %>%
  table(useNA = "always") %>%
  as.data.frame()
motifs_with_SELEX_data %>%
  filter(is.na(type)) %>%
  pull(ID)
motifs_with_SELEX_data$type[which(is.na(motifs_with_SELEX_data$type))] <- motifs_with_SELEX_data$type_review[which(is.na(motifs_with_SELEX_data$type))]

motifs_with_SELEX_data$type %>%
  table(useNA = "always") %>%
  as.data.frame()
#                . Freq
# 1       composite 1439
# 2      composite?    1
# 3         dimeric  808
# 4        dimeric?   22
# 5       monomeric 1051
# 6      monomeric?   11
# 7         spacing  450
# 8        spacing?    8
# 9   ssDNA binding   73
# 10 ssDNA binding?    2
# 11     tetrameric   17
# 12    tetrameric?   10
# 13       trimeric   19
# 14        unknown   22
# 15           <NA>    0

motifs_with_SELEX_data %>%
  pull(type_review) %>%
  table(useNA = "always") %>%
  as.data.frame()

#                . Freq
# 1      composite 1445
# 2        dimeric  828
# 3      monomeric 1052
# 4     monomeric?    2
# 5        spacing  453
# 6  ssDNA binding   73
# 7     tetrameric   21
# 8       trimeric   25
# 9        unknown   22
# 10          <NA>   12


motifs_with_SELEX_data %>%
  filter(is.na(type_review)) %>%
  select(ID, representative)

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  separate(symbol,
    into = c("TF1", "TF2"),
    sep = "_",
    fill = "right",
    remove = FALSE
  ) %>%
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


motifs_with_SELEX_data$sorted_TF2_copies <- NA
# Check that the TF1_sorted and TF2_sorted match between tables

heterodimer_ids <- motifs_with_SELEX_data %>%
  filter(experiment == "CAP-SELEX") %>%
  pull(ID)
row_numbers <- which(motifs_with_SELEX_data$ID %in% heterodimer_ids)

match_inds <- match(heterodimer_ids, heterodimer_nodes$name)

head(heterodimer_ids)
head(heterodimer_nodes$name[match_inds])

table((motifs_with_SELEX_data %>% filter(experiment == "CAP-SELEX") %>% pull(TF1_sorted)) == heterodimer_nodes$TF1_sorted[match_inds]) # All TRUE
table((motifs_with_SELEX_data %>% filter(experiment == "CAP-SELEX") %>% pull(TF2_sorted)) == heterodimer_nodes$TF2_sorted[match_inds]) # All TRUE

motifs_with_SELEX_data$sorted_TF1_copies[row_numbers] <- heterodimer_nodes$sorted_TF1_copies[match_inds]
motifs_with_SELEX_data$sorted_TF2_copies[row_numbers] <- heterodimer_nodes$sorted_TF2_copies[match_inds]

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  relocate(sorted_TF1_copies, sorted_TF2_copies, type_review, .after = TF2_sorted)


# Add TF1_copies and TF2_copies

motifs_with_SELEX_data$TF1_copies <- motifs_with_SELEX_data$sorted_TF1_copies
motifs_with_SELEX_data$TF2_copies <- motifs_with_SELEX_data$sorted_TF2_copies

swap_ids <- row_numbers[which(motifs_with_SELEX_data$TF1[row_numbers] != motifs_with_SELEX_data$TF1_sorted[row_numbers])]

motifs_with_SELEX_data$TF1_copies[swap_ids] <- motifs_with_SELEX_data$sorted_TF2_copies[swap_ids]
motifs_with_SELEX_data$TF2_copies[swap_ids] <- motifs_with_SELEX_data$sorted_TF1_copies[swap_ids]

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  relocate(TF1_copies, TF2_copies, type_review, .after = TF2)

# For how many heterodimeric proteins the number of copies is unknown

motifs_with_SELEX_data %>%
  filter(experiment == "CAP-SELEX" & representative == "YES") %>%
  filter(grepl("\\?", TF1_copies) | grepl("\\?", TF2_copies)) %>%
  nrow() # 254 (129 representative)
motifs_with_SELEX_data %>%
  filter(experiment == "CAP-SELEX") %>%
  nrow() # 1896
motifs_with_SELEX_data %>%
  filter(experiment == "CAP-SELEX" & representative == "YES") %>%
  nrow() # 793
# Are there some cases with more than 4 proteins

table(rowSums(data.frame(X1 = as.numeric(gsub("\\?", "", motifs_with_SELEX_data$TF1_copies)), X2 = as.numeric(gsub("\\?", "", motifs_with_SELEX_data$TF2_copies))), na.rm = TRUE))
#    0    1    2    3    4
#  111 1052 2139  574   57

# Load protein sequence data ----------------------------------------------

data_file <-
  "motifs_with_proteins.tsv"

motifs_with_proteins <- read_delim(
  paste0(
    "motif-metadata/",
    data_file
  ),
  delim = "\t"
)

# Remove protein sequences that are wrong
motifs_with_proteins <- motifs_with_proteins %>%
  mutate(
    TF1_protein_sequence = if_else(TF1_dbd_matches_claimed_gene == FALSE, NA_character_, TF1_protein_sequence),
    TF2_protein_sequence = if_else(TF2_dbd_matches_claimed_gene == FALSE, NA_character_, TF2_protein_sequence)
  )

motifs_with_SELEX_data$`protein sequence 1` <- NA
motifs_with_SELEX_data$`protein sequence 2` <- NA
motifs_with_SELEX_data$`protein sequence 3` <- NA
motifs_with_SELEX_data$`protein sequence 4` <- NA

motifs_with_SELEX_data <- motifs_with_SELEX_data %>% left_join(
  motifs_with_proteins %>%
    select(ID, setdiff(names(motifs_with_proteins), names(motifs_with_SELEX_data))[-1]),
  by = "ID"
)


# Monomers ---------------------------------------------------------------

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  mutate(
    n = case_when(
      type_review == "monomeric" ~ 1,
      type_review == "dimeric" ~ 2,
      type_review == "trimeric" ~ 3,
      type_review == "tetrameric" ~ 4,
      TRUE ~ 0
    ),
    `protein sequence 1` = if_else(n >= 1 & !is.na(TF1_protein_sequence),
      TF1_protein_sequence, `protein sequence 1`
    ),
    `protein sequence 2` = if_else(n >= 2 & !is.na(TF1_protein_sequence),
      TF1_protein_sequence, `protein sequence 2`
    ),
    `protein sequence 3` = if_else(n >= 3 & !is.na(TF1_protein_sequence),
      TF1_protein_sequence, `protein sequence 3`
    ),
    `protein sequence 4` = if_else(n >= 4 & !is.na(TF1_protein_sequence),
      TF1_protein_sequence, `protein sequence 4`
    )
  ) %>%
  select(-n)

# save.image("/scratch/project_2013895/SELEX/RData/tmp.RData")
# load("/scratch/project_2013895/SELEX/RData/tmp.RData")
# Heterodimers ------------------------------------------------------------

## Jolma2015 Heterodimers ------------------------------------------------------------


# Check that the TF1 and TF2 match

match_ind <- match(motifs_with_proteins %>% filter(study == "Jolma2015" & experiment == "CAP-SELEX") %>% pull(ID), motifs_with_SELEX_data$ID)
na_ind <- which(is.na(match_ind))

((motifs_with_proteins %>% filter(study == "Jolma2015" & experiment == "CAP-SELEX") %>% pull(TF1)) == motifs_with_SELEX_data$TF1[match_ind]) %>% table() # all true
((motifs_with_proteins %>% filter(study == "Jolma2015" & experiment == "CAP-SELEX") %>% pull(TF2)) == motifs_with_SELEX_data$TF2[match_ind]) %>% table() # all true


# What kind of cases there can be: 1-1, 1-2, 2-1, 2-2, 1-3

apply(motifs_with_SELEX_data[which(motifs_with_SELEX_data$experiment == "CAP-SELEX"), c("TF1_copies", "TF2_copies")], 1, paste0, collapse = "_") %>% unique()



## Xie2025 Heterodimers ------------------------------------------------------------


# Check that the TF1 and TF2 match

match_ind <- match(motifs_with_proteins %>% filter(study == "Xie2025" & experiment == "CAP-SELEX") %>% pull(ID), motifs_with_SELEX_data$ID)
na_ind <- which(is.na(match_ind))

((motifs_with_proteins %>% filter(study == "Xie2025" & experiment == "CAP-SELEX") %>% pull(TF1)) == motifs_with_SELEX_data$TF1[match_ind]) %>% table() # all true
((motifs_with_proteins %>% filter(study == "Xie2025" & experiment == "CAP-SELEX") %>% pull(TF2)) == motifs_with_SELEX_data$TF2[match_ind]) %>% table() # all true


# What kind of cases there can be: 1-1, 1-2, 2-1, 2-2, 1-3

apply(motifs_with_SELEX_data[which(motifs_with_SELEX_data$experiment == "CAP-SELEX"), c("TF1_copies", "TF2_copies")], 1, paste0, collapse = "_") %>% unique()


motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  rowwise() %>%
  mutate(
    valid = grepl("^\\d+$", TF1_copies) &
      grepl("^\\d+$", TF2_copies) &
      !is.na(TF1_protein_sequence) &
      !is.na(TF2_protein_sequence),
    proteins = list(
      if (valid) {
        x <- c(
          rep(TF1_protein_sequence, as.integer(TF1_copies)),
          rep(TF2_protein_sequence, as.integer(TF2_copies))
        )
        length(x) <- 4
        x
      } else {
        c(
          `protein sequence 1`,
          `protein sequence 2`,
          `protein sequence 3`,
          `protein sequence 4`
        )
      }
    ),
    `protein sequence 1` = proteins[1],
    `protein sequence 2` = proteins[2],
    `protein sequence 3` = proteins[3],
    `protein sequence 4` = proteins[4]
  ) %>%
  ungroup() %>%
  select(
    -proteins, -valid
  )


# Stats -------------------------------------------------------------------

# For how many we were able to find protein sequences

nrow(motifs_with_SELEX_data) # 3933

(is.na(motifs_with_SELEX_data %>% filter(representative == "YES") %>% pull(`protein sequence 1`))) %>% table()
# FALSE  TRUE
# 1046   186

is.na(motifs_with_SELEX_data$`protein sequence 1`) %>% table()
# FALSE  TRUE
# 3247   686

tmp <- motifs_with_SELEX_data %>% select(ID, type_review, TF1, TF2, TF1_copies, TF2_copies, `protein sequence 1`, `protein sequence 2`, `protein sequence 3`, `protein sequence 4`)

# write_delim(motifs_with_SELEX_data,
#            "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv", delim="\t")

# write_delim(motifs_with_SELEX_data,
#             "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25022026.tsv", delim="\t")

# write_delim(motifs_with_SELEX_data,
#   "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_04032026.tsv",
#   delim = "\t"
# )

# Remove some columns that are not needed?

# Add info whether we have the k-mer scores
files <- dir("Data/sequence-to-affinity-data/separate_with_protein_sequences/")
files <- gsub(".tsv", "", files)
motifs_with_SELEX_data$kmer_scores_available <- motifs_with_SELEX_data$motif_ID %in% files

motifs_with_SELEX_data %>%
  pull(kmer_scores_available) %>%
  table()
# FALSE  TRUE
#  826  3107

# add pubmed IDs of the papers

# Jolma2013 PMID: 23332764
# Morgunova2015 PMID: 26632596
# Jolma2015 PMID: 26550823
# Nitta PMID: 25779349
# Xie2025 PMID: 40205063
# Yin2017 PMID: 28473536

motifs_with_SELEX_data <- motifs_with_SELEX_data %>%
  mutate(
    PMID = case_when(
      study == "Jolma2013" ~ "23332764",
      study == "Morgunova2015" ~ "26632596",
      study == "Jolma2015" ~ "26550823",
      study == "Nitta2015" ~ "25779349",
      study == "Xie2025" ~ "40205063",
      study == "Yin2017" ~ "28473536",
      TRUE ~ NA_character_
    )
  )
motifs_with_SELEX_data %>%
  select(PMID) %>%
  unique()

motifs_with_SELEX_data %>%
  filter(!is.na(Human_Ensemble_ID)) %>%
  filter(Human_Ensemble_ID == `TF1 Ensembl_ID`) %>%
  nrow() # 1534

motifs_with_SELEX_data %>%
  filter(!is.na(Human_Ensemble_ID)) %>%
  filter(Human_Ensemble_ID == `TF1_ensembl_gene_id`) %>%
  nrow() # 1509

test <- motifs_with_SELEX_data %>%
  filter(!is.na(Human_Ensemble_ID)) %>%
  filter(Human_Ensemble_ID != `TF1 Ensembl_ID`) %>%
  select(TF1, Human_Ensemble_ID, `TF1 Ensembl_ID`, TF1_alternative_Ensembl_ID, TF1_ensembl_gene_id) # 34

# Remove Methyl-HT-SELEX motifs, ChIP-seq motifs, and ssDNA binding, and unknown motifs

data <- motifs_with_SELEX_data %>% filter(experiment != "Methyl-HT-SELEX" & study != "Morgunova2015" & !(type_review %in% c("ssDNA binding"))) # 3540

data <- data %>% mutate(
  SELEX_data_available = if_else(!is.na(CSC_SELEX_filename) & !is.na(CSC_SELEX_background_filename), "yes", "no")
)

data %>%
  filter(is.na(submitted_ftp_signal)) %>%
  nrow() # 29

data %>%
  filter(is.na(fastq_ftp_signal)) %>%
  nrow() # 38?

data %>%
  filter(is.na(CSC_SELEX_filename)) %>%
  nrow() # 29

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
# 38 3524

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
# 27 3535

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
# 288 3274
# stoichiometry_available,kmer_scores_available,protein_sequences_available,SELEX_data_available

# Draw Venn diagram of four sets: SELEX data available, k-mer scores available, protein sequences available, stoichiometry available

library(ggVennDiagram)
library(UpSetR)

venn_sets <- list(
  "SELEX data\navailable" = data$ID[data$SELEX_data_available == "yes"],
  "k-mer scores\navailable" = data$ID[data$kmer_scores_available == TRUE],
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
    `k-mer scores available` = as.integer(kmer_scores_available),
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

data <- data %>% filter(kmer_scores_available == TRUE & SELEX_data_available == "yes" & protein_sequences_available == "yes" & stoichiometry_available == "yes") # 3095

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


# What columns to keep in the table sent to the community?

cols_to_keep <- c(
  "motif_ID",
  "TF1", "TF1_ensembl_gene_id", "organism1", "TF1_copies", "TF1_protein_sequence",
  "TF2", "TF2_ensembl_gene_id", "organism2", "TF2_copies", "TF2_protein_sequence",
  "protein sequence 1", "protein sequence 2", "protein sequence 3", "protein sequence 4",
  "type_review",
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
  "kmer_scores_available", "SELEX_data_available",
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
  select(type_review) %>%
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

data <- data %>% select(-TF1_dbd_matches_claimed_gene, -TF2_dbd_matches_claimed_gene, -kmer_scores_available, -SELEX_data_available, -protein_sequences_available, -stoichiometry_available)

sum(data$read_count_signal, na.rm = TRUE) + sum(data$read_count_background, na.rm = TRUE) # 1 720 602 884 ~ 1.7 billion reads

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

write_delim(data,
  "motif-metadata/S2A_curated_01042026.tsv",
  delim = "\t"
)

write_delim(motifs_with_SELEX_data,
  "motif-metadata/motifs_with_SELEX_signal_and_background_and_proteins_01042026.tsv",
  delim = "\t"
)
