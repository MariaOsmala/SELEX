library("readr")
library("tidyverse")
library("readxl")

rm(list = ls())
df_motif_info <- read_delim(
  "~/projects/TFBS/Data/SELEX-motif-collection/metadata_final.tsv", "\t"
)

# Jolma2013 ---------------------------------------------------------------

df_Jolma2013 <- df_motif_info %>% filter(study == "Jolma2013")

df_Jolma2013=df_Jolma2013 %>% rename(TF1=symbol)

# Replace instances of FL in column clone by "full"

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`clone` = ifelse(`clone` == "FL", "full", `clone`))

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`TF1` = ifelse(
    `clone` == "mouse DBD_mutant_DBD",
    "Egr1_E410D_FARSDERtoFARSDDR", `TF1`
  ))

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(`clone` = ifelse(`clone` == "mouse DBD_mutant_DBD", "DBD", `clone`))

df_Jolma2013 <- df_Jolma2013 %>%
  mutate(TF1_clone = paste0(TF1, "_", clone)) %>%
  relocate(TF1_clone, .after = clone)


## Read protein sequence info ----------------------------------------------

# Jolma 2013: Table S1 Protein and DNA sequence information, related to Figure 1

# A: Sequences of proteins used and their source;
# the amino-acid sequence is based on end-sequencing of DNA and ENSEMBL.
#
# Clones are either
#   synthesized (Gene synthesis),
#   PCR-cloned from the indicated source,
#   or Gateway-cloned from the human ORFeome (hORFeome).

protein_sequences <- read_excel(
  "~/projects/SELEX/TF-amino-acid-sequences/Data/1-s2.0-S0092867412014961-mmc1.xls",
  skip = 3,
  n_max = 539
)

# HNGC-name"            "Ensembl id."
# "amino-acid sequence" "main structural class" "source organism" "clone type (DBD, full)"
# "production method"  "Clone source"

protein_sequences$`TF1 Classificaton in Vaquerizas et al., 2009` <- NA

# Try to equalise the names

protein_sequences <- protein_sequences %>%
  rename(`TF1 HNGC` = `HNGC-name`)

protein_sequences <- protein_sequences %>%
  rename(`TF1 Ensembl_ID` = `Ensembl id.`)

protein_sequences <- protein_sequences %>%
  rename(`TF1_protein_sequence` = `amino-acid sequence`)

protein_sequences <- protein_sequences %>%
  rename(`TF1 DOMAIN` = `main structural class`,
  "TF1 source organism" ="source organism",
   "TF1 clone type" = "clone type",
  "TF1 production method"="production method",
  "TF1 Clone source"="Clone source")

# Correct some protein names
protein_sequences <- protein_sequences %>%
  mutate(`TF1 HNGC` = ifelse(`TF1 HNGC` == "Trp53", "Tp53", `TF1 HNGC`))

protein_sequences <- protein_sequences %>%
  mutate(`TF1 HNGC` = ifelse(`TF1 HNGC` == "Trp73", "Tp73", `TF1 HNGC`))

protein_sequences <- protein_sequences %>%
  mutate(TF1_clone = paste0(`TF1 HNGC`, "_", `TF1 clone type`))

# Are the protein sequences unique for the same symbol_clone?
protein_sequences$TF1_clone %>% length() # 538
protein_sequences$TF1_clone %>%
  unique() %>%
  length() # 538
protein_sequences$`TF1_protein_sequence` %>%
  unique() %>%
  length() # 534

df_Jolma2013_protein <- df_Jolma2013 %>%
  left_join(protein_sequences, by = "TF1_clone")

# Did we find protein for all motifs YES!
df_Jolma2013_protein$`TF1_protein_sequence` %>%
  is.na() %>%
  table(useNA = "always")

# Ensembl IDS no not always match and some missing
(df_Jolma2013_protein$Human_Ensemble_ID == df_Jolma2013_protein$`TF1 Ensembl_ID`) %>% table(useNA = "always")
# FALSE  TRUE  <NA>
#  20   654   146

# 146 missing because it was not in Lambert table, are these mouse motifs. YES and some human motifs.
df_Jolma2013_protein$Human_Ensemble_ID %>%
  is.na() %>%
  table()
# FALSE  TRUE
#  674   146

df_Jolma2013_protein %>%
  filter(is.na(Human_Ensemble_ID) & organism=="Homo_sapiens") %>%
  pull(TF1) %>% unique()
# "HINFP1"   "ZNF238"   "ZNF306"   "POU5F1P1" "BHLHB2"   "BHLHB3"   "CART1"    "RAXL1"    "ZNF435"

df_Jolma2013_protein %>%
  filter(is.na(Human_Ensemble_ID)) %>%
  select(organism) %>%
  table()
# Homo_sapiens Mus_musculus
# 13          133

names(df_Jolma2013_protein)


Jolma2013_protein_sequences_unique <- protein_sequences

# Try to add info to protein sequences (538 x 9)
Jolma2013_protein_sequences <- protein_sequences %>% # 820 x 25
  left_join(df_Jolma2013 %>%
    select(
      ID, TF1_clone, study, experiment, family, Lambert2018_families,
      ligand, batch, seed, multinomial, cycle, representative, type, comment,
      filename, Methyl.SELEX.Motif.Category, Human_Ensemble_ID
    ), by = "TF1_clone")





# Yin 2017 --------------------------------------------------------------

df_Yin2017 <- df_motif_info %>% filter(study == "Yin2017" & experiment == "HT-SELEX")
df_Yin2017_Methyl <- df_motif_info %>% filter(study == "Yin2017" & experiment == "Methyl-HT-SELEX")

df_Yin2017=df_Yin2017 %>% rename(TF1=symbol)
df_Yin2017_Methyl=df_Yin2017_Methyl %>% rename(TF1=symbol)


df_Yin2017 <- df_Yin2017 %>%
  mutate(TF1_clone = paste0(TF1, "_", clone)) %>%
  relocate(TF1_clone, .after = clone)

df_Yin2017_Methyl <- df_Yin2017_Methyl %>%
  mutate(TF1_clone = paste0(TF1, "_", clone)) %>%
  relocate(TF1_clone, .after = clone)


# Table S1 - Sequence information for Proteins (Section A, top)

# Section a
# This table contains information about all of the TF constructs used in current study.
# Columns from left to right contain:
# Ensembl id.	Gene identifier
# HNGC-name	HUGO Gene Nomenclature Committee accepted symbol for the TF
# Construct type	Full length (FL) or extended DNA Binding Domain (eDBD)
# main structural class	Main DNA binding domain in the construct
# amino-acid sequence
# Classification in Ref. Vaquerizas et al., 2009


## Read protein sequence info ----------------------------------------------

protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/aaj2239_yin_sm_tables_s1-s6.xlsx",
  sheet = "S1 sequence information", skip = 17, n_max = 1178 - 17
)

names(protein_sequences)

# Try to equalise the names

protein_sequences <- protein_sequences %>%
  rename(`TF1 HNGC` = `HNGC`)

protein_sequences <- protein_sequences %>%
  rename(`TF1 Ensembl_ID` = `Ensembl_ID`)

protein_sequences$`TF1 source organism` <- "human"
protein_sequences$`TF1 Clone source`=NA
protein_sequences$`TF1 production method`=NA

protein_sequences <- protein_sequences %>%
  rename(`TF1 DOMAIN` = `DOMAIN`,
  `TF1 clone type` = `Clone used`,
  `TF1 Classificaton in Vaquerizas et al., 2009` = `Classificaton in Vaquerizas et al., 2009`
    )


doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))
#"GLIS2"  "HOXA10" "MEF2B"  "PBX2"   "SPIB"
ind=which(protein_sequences$`TF1 HNGC`==doubles[1])
protein_sequences$`eDBD seuqence`[ind[1]]==protein_sequences$`eDBD seuqence`[ind[2]] #TRUE

protein_sequences$TF1_alternative_Ensembl_ID=NA
protein_sequences$TF1_alternative_Ensembl_ID[ind[1]]=protein_sequences$`TF1 Ensembl_ID`[ind[2]]
protein_sequences=protein_sequences[-ind[2],]


doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))
ind=which(protein_sequences$`TF1 HNGC`==doubles[1])
test=protein_sequences[ind, ]

protein_sequences$`TF1 clone type`[ind[1]]="eDBD; FL"
protein_sequences$`Full-length sequence`[ind[1]]=protein_sequences$`Full-length sequence`[ind[2]]
protein_sequences$TF1_alternative_Ensembl_ID[ind[1]]=protein_sequences$`TF1 Ensembl_ID`[ind[2]]
protein_sequences=protein_sequences[-ind[2],]

doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))

ind=which(protein_sequences$`TF1 HNGC`==doubles[1])
protein_sequences$`eDBD seuqence`[ind[1]]==protein_sequences$`eDBD seuqence`[ind[2]] #same sequence
#" RKKIQISRILDQRNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSANRLFQYASTDMDRVLLKYTEYSEPHESRTNTDILEVPQTLKRRGIGLDGPELEPDE"
#"GRKKIQISRILDQRNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSANRLFQYASTDMDRVLLKYTEYSEPHESRTNTDILEVPQTLKRRGIGLDGPELEPDE"
protein_sequences$TF1_alternative_Ensembl_ID[ind[2]]=protein_sequences$`TF1 Ensembl_ID`[ind[1]]
protein_sequences=protein_sequences[-ind[1],]


doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))

ind=which(protein_sequences$`TF1 HNGC`==doubles[1])
protein_sequences$`eDBD seuqence`[ind[1]]==protein_sequences$`eDBD seuqence`[ind[2]] #same sequence
#" RKKIQISRILDQRNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSANRLFQYASTDMDRVLLKYTEYSEPHESRTNTDILEVPQTLKRRGIGLDGPELEPDE"
#"GRKKIQISRILDQRNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSANRLFQYASTDMDRVLLKYTEYSEPHESRTNTDILEVPQTLKRRGIGLDGPELEPDE"
protein_sequences$TF1_alternative_Ensembl_ID[ind[2]]=protein_sequences$`TF1 Ensembl_ID`[ind[1]]
protein_sequences=protein_sequences[-ind[1],]


doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))

ind=which(protein_sequences$`TF1 HNGC`==doubles[1])
test=protein_sequences[ind, ]
protein_sequences$`eDBD seuqence`[ind[1]]==protein_sequences$`eDBD seuqence`[ind[2]] #same sequence
protein_sequences$TF1_alternative_Ensembl_ID[ind[2]]=protein_sequences$`TF1 Ensembl_ID`[ind[1]]
protein_sequences$`TF1 clone type`[ind[2]]="eDBD; FL"
protein_sequences=protein_sequences[-ind[1],]

doubles=names(which(protein_sequences$`TF1 HNGC` %>% table()>1))

protein_sequences %>% select(`TF1 Classificaton in Vaquerizas et al., 2009`) %>% table(useNA = "always")
# Class a
# Class a, ensembl ID is different
# Class b
# Class b, ensembl ID is different
# Class c
# Class c, ensembl ID is different
# Class other
# Class x
# Class x, ensembl ID is different
# Not in Vaquerizas et al list
# <NA>

splits =strsplit(protein_sequences %>% pull(`TF1 Classificaton in Vaquerizas et al., 2009`), ", ")
lengths <- sapply(splits, length)
table(lengths)

max_len <- max(sapply(splits, length))

# Pad with NAs
padded <- lapply(splits, function(x) {
  length(x) <- max_len
  return(x)
})

comment <- as.data.frame(do.call(rbind, padded), stringsAsFactors = FALSE)

protein_sequences$`TF1 Classificaton in Vaquerizas et al., 2009` <- comment$V1
protein_sequences$`TF1 Vaquerizas_comment` <- comment$V2

protein_sequences <- protein_sequences %>%
  # Fix the column name typo
  rename(`TF1 eDBD sequence` = `eDBD seuqence`,
  `TF1 Full-length sequence` = `Full-length sequence`
  ) %>%
  separate_rows(`TF1 clone type`, sep = ";\\s*") # split into separate rows on ";"

protein_sequences <- protein_sequences %>%
  mutate(
    `TF1 eDBD sequence` = ifelse(`TF1 clone type` == "eDBD", `TF1 eDBD sequence`, NA),
    `TF1 Full-length sequence` = ifelse(`TF1 clone type` == "FL", `TF1 Full-length sequence`, NA)
  ) %>%
  mutate(TF1_protein_sequence = coalesce(`TF1 eDBD sequence`, `TF1 Full-length sequence`)) %>%
  # Optionally remove the original two columns
  select(-`TF1 eDBD sequence`, -`TF1 Full-length sequence`)

protein_sequences <- protein_sequences %>%
  mutate(TF1_clone = paste0(`TF1 HNGC`, "_", `TF1 clone type`))


names(protein_sequences)

# Test for whitespaces in the symbol_clone column

protein_sequences$TF1_clone <- str_replace_all(protein_sequences$TF1_clone, "\\s+", "")

# There is * is protein sequences? Remove it
protein_sequences$TF1_protein_sequence[grep("\\*", protein_sequences$TF1_protein_sequence)]

protein_sequences$TF1_protein_sequence <- gsub("\\*", "", protein_sequences$TF1_protein_sequence)

# which(str_detect(df_Yin2017$TF1_clone, "\\s") )
# which(str_detect(protein_sequences$TF1_clone, "\\s") )


df_Yin2017_protein <- df_Yin2017 %>%
  left_join(protein_sequences, by = "TF1_clone")

test <- df_Yin2017_protein[df_Yin2017_protein$TF1_protein_sequence %>%
  is.na() %>%
  which(), ]

test %>% select(TF1_clone, Human_Ensemble_ID, `TF1 Ensembl_ID`)

# TF1_clone Human_Ensemble_ID Ensembl_ID
#<chr>        <chr>             <chr>
# 1 SIX2_eDBD    ENSG00000170577   ENSG00000170577                           No eDBD protein, only full

# Maybe eDBD and full are the same?

# The following have different names in the protein_sequences table,
# the ensembl IDs are the same

# 3 FOXO3_eDBD   ENSG00000118689   NA                 ENSG00000118689 FOXO3A
# 4 FOXG1_eDBD   ENSG00000176165   NA                 ENSG00000176165 FOXG1B
# 6 SKOR1_eDBD   ENSG00000188779   NA                 ENSG00000188779 LBXCOR1
# 7 ZSCAN5A_eDBD ENSG00000131848   NA                 ENSG00000131848 ZSCAN5


missing_ind <- df_Yin2017_protein$TF1_protein_sequence %>%
  is.na() %>%
  which()

df_Yin2017_protein$TF1[missing_ind]

lookup <- c(FOXO3 = "FOXO3A", FOXG1 = "FOXG1B", SKOR1 = "LBXCOR1", ZSCAN5A = "ZSCAN5")

for (i in 1:length(lookup)) {
  missing_ind <- which(df_Yin2017_protein$TF1 == names(lookup[i]))

  # test=df_Yin2017_protein[missing_ind, ]

  df_Yin2017_protein$TF1_clone[missing_ind] <- gsub(names(lookup[i]), lookup[i], df_Yin2017_protein$TF1_clone[missing_ind])
  df_Yin2017$TF1_clone[missing_ind] <- gsub(names(lookup[i]), lookup[i], df_Yin2017$TF1_clone[missing_ind])

  df_Yin2017_protein[missing_ind, ] <- df_Yin2017_protein[missing_ind, 1:ncol(df_Yin2017)] %>%
    left_join(protein_sequences, by = "TF1_clone")
}

test <- df_Yin2017_protein[df_Yin2017_protein$TF1_protein_sequence %>%
  is.na() %>%
  which(), ]

test %>% select(TF1_clone, Human_Ensemble_ID, `TF1 Ensembl_ID`)
# 1 SIX2_eDBD    ENSG00000170577   ENSG00000170577
names(df_Yin2017_protein)
# "TF1 Ensembl_ID"                               "TF1 HNGC"                                     "TF1 clone type"
# "TF1 DOMAIN"                                   "TF1 Classificaton in Vaquerizas et al., 2009" "TF1 source_organism"
# "TF1 Clone source"                             "TF1 production method"                        "TF1_alternative_Ensembl_ID"
# "TF1 Vaquerizas_comment"                       "TF1_protein_sequence"

#Does df_Jolma2013 and df_Yin2017 have the same columns

setdiff(names(df_Jolma2013_protein), names(df_Yin2017_protein)) #empty
setdiff( names(df_Yin2017_protein), names(df_Jolma2013_protein)) #empty

df_Jolma2013_protein$TF1_alternative_Ensembl_ID=NA
df_Jolma2013_protein$`TF1 Vaquerizas_comment`=NA

# Is there some for which the protein info is missing:

df_Yin2017_protein$TF1_protein_sequence %>%
  is.na() %>%
  table(useNA = "always")
# FALSE  TRUE  <NA>
#  863     1    0

df_Yin2017_protein %>%
  filter(is.na(TF1_protein_sequence)) %>%
  select(ID, TF1_clone)
# SIX2_HT-SELEX_TTACTA40NAGG_KV_NCGTATCRYN_1_4 SIX2_eDBD

# Try to add info to protein sequences (1534    7)

Yin2017_protein_sequences_unique <- protein_sequences # 1530

Yin2017_protein_sequences <- protein_sequences %>% # 1724
  left_join(df_Yin2017 %>%
    select(
      ID, TF1_clone, study, experiment, family, Lambert2018_families,
      ligand, batch, seed, multinomial, cycle, representative, type, comment,
      filename, Methyl.SELEX.Motif.Category, Human_Ensemble_ID
    ), by = "TF1_clone")


# How many motifs there are for each protein
# Yin2017_protein_sequences$number_of_motifs=sapply(Yin2017_protein_sequences$`df_Yin2017 %>% ...`, nrow)

Yin2017_protein_sequences %>%
  filter(!is.na(ID)) %>%
  nrow() # 864

Yin2017_protein_sequences %>%
  filter(!is.na(ID)) %>%
  pull(TF1_clone) %>%
  unique() %>%
  length() # 670

tmp <- Yin2017_protein_sequences %>%
  filter(!is.na(ID)) %>%
  select(`TF1 HNGC`, `TF1 clone type`)

paste0(tmp$`TF1 HNGC`, "_", tmp$`TF1 clone type`) %>%
  unique() %>%
  length() # 670


setdiff(names(Jolma2013_protein_sequences), names(Yin2017_protein_sequences))
setdiff(names(Yin2017_protein_sequences), names(Jolma2013_protein_sequences))

Jolma2013_protein_sequences$TF1_clone %in% gsub("eDBD", "DBD", Yin2017_protein_sequences$TF1_clone) %>% table()
# FALSE  TRUE
# 375   445

gsub("eDBD", "DBD", Yin2017_protein_sequences$TF1_clone) %in% Jolma2013_protein_sequences$TF1_clone %>% table()
# FALSE  TRUE
# 1357   367


# Are the protein sequences same in both studies

match_ind <- match(Jolma2013_protein_sequences$TF1_clone, gsub("eDBD", "DBD", Yin2017_protein_sequences$TF1_clone))
na.ind <- is.na(match_ind)

head(Jolma2013_protein_sequences$TF1_clone[-which(na.ind)], 10)
head(gsub("eDBD", "DBD", Yin2017_protein_sequences$TF1_clone)[match_ind[-which(na.ind)]], 10)

(Jolma2013_protein_sequences$TF1_protein_sequence[-which(na.ind)] == Yin2017_protein_sequences$TF1_protein_sequence[match_ind[-which(na.ind)]]) %>% table(useNA = "always")
# FALSE  TRUE  <NA>
# 153   292   0

# Some are, some are not, maybe changes are not big ones

not_equal <- which(!(Jolma2013_protein_sequences$TF1_protein_sequence[-which(na.ind)] == Yin2017_protein_sequences$TF1_protein_sequence[match_ind[-which(na.ind)]]))

Jolma2013_protein_sequences[which(!na.ind)[not_equal[1]], c("TF1_clone", "TF1_protein_sequence")]
Yin2017_protein_sequences[match_ind[which(!na.ind)][not_equal[1]], c("TF1_clone", "TF1_protein_sequence")]

## Add Vaquerizas classification to Jolma2013 motifs

Vaquerizas_info <- Yin2017_protein_sequences %>%
  dplyr::group_by(`TF1 HNGC`) %>%
  dplyr::summarise(Va = unique(`TF1 Classificaton in Vaquerizas et al., 2009`), .groups = "drop")

match_ind <- match(Jolma2013_protein_sequences$`TF1 HNGC`, Vaquerizas_info$`TF1 HNGC`)
na.ind <- is.na(match_ind)

head(Jolma2013_protein_sequences$`TF1 HNGC`[-which(na.ind)], 50)
head(Vaquerizas_info$`TF1 HNGC`[match_ind[-which(na.ind)]], 50)

Jolma2013_protein_sequences$`TF1 Classificaton in Vaquerizas et al., 2009`[-which(na.ind)] <- Vaquerizas_info$Va[match_ind[-which(na.ind)]]


#
Jolma2013_protein_sequences$ID %>%
  unique() %>%
  length() # 820

nrow(Yin2017_protein_sequences) # 1724

Yin2017_protein_sequences$ID %>%
  unique() %>%
  length() # 865
Yin2017_protein_sequences %>%
  filter(!is.na(ID)) %>%
  pull(ID) %>%
  unique() %>%
  length() # 864

# Add Vaqueriza classification to Jolma2013 motifs
df_Jolma2013_protein$TF1_protein_sequence %>%
  is.na() %>%
  table(useNA = "always") # 820 of 820

df_Yin2017_protein$TF1_protein_sequence %>%
  is.na() %>%
  table(useNA = "always") # 863 has protein, 1 missing

df_Jolma2013_protein %>%
  filter(TF1 == "E2F2") %>%
  select(TF1_clone, TF1_protein_sequence)
df_Yin2017_protein %>%
  filter(TF1 == "E2F2") %>%
  select(TF1_clone, TF1_protein_sequence)


# Methyl SELEX motifs --------------------------------------------------------

df_Yin2017_Methyl_protein <- df_Yin2017_Methyl %>%
  left_join(protein_sequences, by = "TF1_clone")

test <- df_Yin2017_Methyl_protein[df_Yin2017_Methyl_protein$TF1_protein_sequence %>%
  is.na() %>%
  which(), ]

test %>% select(TF1_clone, Human_Ensemble_ID, `TF1 Ensembl_ID`)


# 4 FOXG1_eDBD   ENSG00000176165   NA                 ENSG00000176165 FOXG1B


missing_ind <- df_Yin2017_Methyl_protein$TF1_protein_sequence %>%
  is.na() %>%
  which()

df_Yin2017_Methyl_protein$TF1[missing_ind]

lookup <- c(FOXG1 = "FOXG1B")

for (i in 1:length(lookup)) {
  missing_ind <- which(df_Yin2017_Methyl_protein$TF1 == names(lookup[i]))

  # test=df_Yin2017_protein[missing_ind, ]

  df_Yin2017_Methyl_protein$TF1_clone[missing_ind] <- gsub(names(lookup[i]), lookup[i], df_Yin2017_Methyl_protein$TF1_clone[missing_ind])
  df_Yin2017_Methyl$TF1_clone[missing_ind] <- gsub(names(lookup[i]), lookup[i], df_Yin2017_Methyl$TF1_clone[missing_ind])

  df_Yin2017_Methyl_protein[missing_ind, ] <- df_Yin2017_Methyl_protein[missing_ind, 1:ncol(df_Yin2017_Methyl)] %>%
    left_join(protein_sequences, by = "TF1_clone", relationship = "many-to-many")
}

test <- df_Yin2017_Methyl_protein[df_Yin2017_Methyl_protein$TF1_protein_sequence %>%
  is.na() %>%
  which(), ]

test %>% select(TF1_clone, Human_Ensemble_ID, `TF1 Ensembl_ID`)

names(df_Yin2017_Methyl_protein)
# "TF1 Ensembl_ID"                               "TF1 HNGC"                                     "TF1 clone type"
# "TF1 DOMAIN"                                   "TF1 Classificaton in Vaquerizas et al., 2009" "TF1 source organism"
# "TF1 Clone source"                             "TF1 production method"                        "TF1_alternative_Ensembl_ID"
# "TF1 Vaquerizas_comment"                       "TF1_protein_sequence"



# Is there some for which the protein info is missing:

df_Yin2017_Methyl_protein$TF1_protein_sequence %>%
  is.na() %>%
  table(useNA = "always") # NO


# Nitta2015 motifs --------------------------------------------------------

#The paper does not report the protein sequences

# Are the proteins in earlier studies

df_Nitta2015 <- df_motif_info %>% filter(study == "Nitta2015")

df_Nitta2015$symbol <- str_replace_all(df_Nitta2015$symbol, "\\s+", "")

# Clone is missing

# Human HT-SELEX data were generated for PAX3, PGR, NR1I2, and NR1I3 using E. coli expressed DBDs,

df_Nitta2015$symbol
# "SREBF2"  "PAX4"    "PAX4"    "PAX3"    "PAX3"    "PGR"(not)     "ONECUT2" "ONECUT1" "NR1I2"   "NR1I3"

without_clone <- c("SREBF2", "PAX4", "ONECUT2", "ONECUT1")

df_Nitta2015$symbol %in% Jolma2013_protein_sequences$`TF1 HNGC`
df_Nitta2015$symbol %in% Yin2017_protein_sequences$`TF1 HNGC`

Jolma2013_protein_sequences_unique %>%
  filter(`TF1 HNGC` %in% without_clone) %>%
  select(`TF1 HNGC`, `TF1_protein_sequence`, "TF1 clone type")
Yin2017_protein_sequences_unique %>%
  filter(`TF1 HNGC` %in% without_clone) %>%
  select(`TF1 HNGC`, `TF1_protein_sequence`, "TF1 clone type")

"ENSG00000082175" %in% Jolma2013_protein_sequences$`TF1 Ensembl_ID`
"ENSG00000082175" %in% Yin2017_protein_sequences$`TF1 Ensembl_ID`

"ENSG00000082175" %in% Yin2017_protein_sequences$TF1_alternative_Ensembl_ID

df_Nitta2015 <- df_Nitta2015 %>% rename(TF1=symbol)

# Jolma2015 ---------------------------------------------------------------

# Clone info missing?

# symbol(s)	HGNC symbol(s) (gene names);
# Either for reference models that were generated for individual TFs in this study or

# for the co-operative pairs, in which case the TF names are separated with an underscore "_".
# First of the TFs is TF1 that was expressed as a SBP tagged clone and was used in first of the affinity separations
# and the second one is the 3xFLAG tagged TF2

df_Jolma2015_monomers <- df_motif_info %>% filter(study == "Jolma2015", experiment == "HT-SELEX") # 31
df_Jolma2015_heterodimers <- df_motif_info %>% filter(study == "Jolma2015", experiment == "CAP-SELEX") # 562

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

protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/41586_2015_BFnature15518_MOESM33_ESM.xlsx",
  sheet = "S1 sequence information", skip = 24, n_max = 237 - 24
)

names(protein_sequences)

protein_sequences$`Construct type` %>% table()
# Crystal Full length protein used in HT-SELEX with mixed proteins
# 2                                                       11
# TF1(SBP)                                              TF2(3xFLAG)
# 105                                                       94

# Remove Crustal and "Full length protein used in HT-SELEX with mixed proteins"

protein_sequences <- protein_sequences %>%
  filter(`Construct type` %in% c("TF1(SBP)", "TF2(3xFLAG)"))


# Are the protein sequences equal for replicate TFs (TF1(SBP)  and TF2(3xFLAG))?
multiple <- names(which(protein_sequences$`HNGC-name` %>% table() > 1))

# The sequences are the same are differ by 1 or 2 amino acids, consider the shorter

for (mtp in multiple) {
  # mtp=multiple[5]
  # test=protein_sequences %>% filter(`HNGC-name`==mtp)
  seq <- protein_sequences %>%
    filter(`HNGC-name` == mtp) %>%
    pull(`amino-acid sequence`)
  x <- nchar(seq)
  # print(abs(outer(x, x, "-")))
  y <- combn(x, 2, function(v) abs(diff(v)))
  x_no_zeros <- y[y != 0]
  if (is.na(mean(x_no_zeros))) {
    # no difference in length,  does the sequences differ, No they are all sam
    print(mtp)
    # print(all(seq ==seq[1]))
    # take the first
    tmp <- protein_sequences %>% filter(`HNGC-name` == mtp)
    tmp <- tmp[1, ]
    protein_sequences <- protein_sequences %>% filter(`HNGC-name` != mtp)
    protein_sequences <- rbind(protein_sequences, tmp)
  } else {
    # print(mtp)
    # print(mean(x_no_zeros)) #if the mean is 0, then all sequences are the same length
    # select some

    tmp <- protein_sequences %>% filter(`HNGC-name` == mtp)
    tmp <- tmp[which.min(nchar(tmp$`amino-acid sequence`)), ]

    protein_sequences <- protein_sequences %>% filter(`HNGC-name` != mtp)
    protein_sequences <- rbind(protein_sequences, tmp)
  }
  # print(mean(x_no_zeros))
}


protein_sequences$symbol <- protein_sequences$`HNGC-name`


df_Jolma2015_monomers_protein <- df_Jolma2015_monomers %>%
  left_join(protein_sequences, by = "symbol")

names(df_Jolma2015_monomers_protein)



df_Jolma2015_monomers_protein <- df_Jolma2015_monomers_protein %>%
  rename(
    `TF1 HNGC` = `HNGC-name`,
    `TF1 Ensembl_ID` = `Ensembl id.`,
    `TF1_protein_sequence` = `amino-acid sequence`,
    `TF1 DOMAIN` = `main structural class`,
    `TF1 6-mer` = `6-mer`,
    `TF1 source organism` = `source organism`,
    `TF1 Clone source` = `Clone source`,
    `TF1 Construct type` = `Construct type`,
    `TF1 Protein amount` = `Protein amount`,
    `TF1 HT-SELEX validation` = `HT-SELEX validation`,
    `TF1 Activity in CAP-SELEX` = `Activity in CAP-SELEX`,
    `TF1 Included domain Ids` = `Included domain Ids`,
    `TF1 Included domain description` = `Included domain description`,
    `TF1 Excluded domain Ids` = `Excluded domain Ids`,
    `TF1 Excluded domain description` = `Excluded domain description`
  )




df_Jolma2015_heterodimers <- df_Jolma2015_heterodimers %>%
  rename("symbol_both" = "symbol")

df_Jolma2015_heterodimers <- df_Jolma2015_heterodimers %>%
  separate(symbol_both, into = c("symbol", "TF2"), sep = "_", remove = FALSE)

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers %>%
  left_join(protein_sequences, by = "symbol")

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename(
    `TF1 HNGC` = `HNGC-name`,
    `TF1 Ensembl_ID` = `Ensembl id.`,
    `TF1_protein_sequence` = `amino-acid sequence`,
    `TF1 DOMAIN` = `main structural class`,
    `TF1 6-mer` = `6-mer`,
    `TF1 source organism` = `source organism`,
    `TF1 Clone source` = `Clone source`,
    `TF1 Construct type` = `Construct type`,
    `TF1 Protein amount` = `Protein amount`,
    `TF1 HT-SELEX validation` = `HT-SELEX validation`,
    `TF1 Activity in CAP-SELEX` = `Activity in CAP-SELEX`,
    `TF1 Included domain Ids` = `Included domain Ids`,
    `TF1 Included domain description` = `Included domain description`,
    `TF1 Excluded domain Ids` = `Excluded domain Ids`,
    `TF1 Excluded domain description` = `Excluded domain description`
  )

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename("TF1" = "symbol")

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename("symbol" = "TF2")

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  left_join(protein_sequences, by = "symbol")

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename(
    `TF2 HNGC` = `HNGC-name`,
    `TF2 Ensembl_ID` = `Ensembl id.`,
    `TF2_protein_sequence` = `amino-acid sequence`,
    `TF2 DOMAIN` = `main structural class`,
    `TF2 6-mer` = `6-mer`,
    `TF2 source organism` = `source organism`,
    `TF2 Clone source` = `Clone source`,
    `TF2 Construct type` = `Construct type`,
    `TF2 Protein amount` = `Protein amount`,
    `TF2 HT-SELEX validation` = `HT-SELEX validation`,
    `TF2 Activity in CAP-SELEX` = `Activity in CAP-SELEX`,
    `TF2 Included domain Ids` = `Included domain Ids`,
    `TF2 Included domain description` = `Included domain description`,
    `TF2 Excluded domain Ids` = `Excluded domain Ids`,
    `TF2 Excluded domain description` = `Excluded domain description`
  )

df_Jolma2015_heterodimers_protein <- df_Jolma2015_heterodimers_protein %>%
  rename("TF2" = "symbol")

# Everything ok

Jolma2015_protein_sequences <- protein_sequences

# Xie 2025 ----------------------------------------------------------------

df_Xie2025_monomers <- df_motif_info %>% filter(study == "Xie2025" & experiment == "HT-SELEX") # 12

# Is HHEX the same as THHEX?
df_Xie2025_monomers[which(df_Xie2025_monomers$symbol == "THHEX"), "symbol"] <- "HHEX"



df_Xie2025_heterodimers <- df_motif_info %>% filter(study == "Xie2025" & experiment == "CAP-SELEX") # 1336

# These are unique proteins
protein_sequences <- read_excel("~/projects/SELEX/TF-amino-acid-sequences/Data/Xie et al. 2025 Supplementary_Tables.xlsx",
  sheet = "Table S1TF and ligand sequences",
  skip = 1270
)
names(protein_sequences)

protein_sequences$`Clone used` <- gsub("eDBD;", "eDBD", protein_sequences$`Clone used`)

protein_sequences$`Clone used` %>% table()
# eDBD   FL
# 430    3

protein_sequences$symbol <- protein_sequences$HNGC

df_Xie2025_monomers_protein <- df_Xie2025_monomers %>%
  left_join(protein_sequences, by = "symbol")

df_Xie2025_monomers_protein <- df_Xie2025_monomers_protein %>%
  rename(
    `TF1 HNGC` = `HNGC`,
    `TF1 Ensembl_ID` = `Ensembl_ID`,
    `TF1_protein_sequence` = `eDBD seuqence`,
    `TF1 DOMAIN` = `DOMAIN`,
    `TF1 clone type` = `Clone used`
      )



df_Xie2025_heterodimers <- df_Xie2025_heterodimers %>%
  rename("symbol_both" = "symbol")

df_Xie2025_heterodimers <- df_Xie2025_heterodimers %>%
  separate(symbol_both, into = c("symbol", "TF2"), sep = "_", remove = FALSE)

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers %>%
  left_join(protein_sequences, by = "symbol")

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename(
    `TF1 Ensembl_ID` = `Ensembl_ID`,
    `TF1 HNGC` = HNGC,
    `TF1 clone type` = `Clone used`,
    `TF1 DOMAIN` = DOMAIN,
    `TF1_protein_sequence` = `eDBD seuqence`
  )

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename("TF1" = "symbol")

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename("symbol" = "TF2")

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  left_join(protein_sequences, by = "symbol")

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename(
    `TF2 Ensembl_ID` = `Ensembl_ID`,
    `TF2 HNGC` = HNGC,
    `TF2 clone type` = `Clone used`,
    `TF2 DOMAIN` = DOMAIN,
    `TF2_protein_sequence` = `eDBD seuqence`
  )

df_Xie2025_heterodimers_protein$`TF2 HNGC` %>% table(useNA = "always") # 218
df_Xie2025_heterodimers_protein$`TF1 HNGC` %>% table(useNA = "always") # 218

df_Xie2025_heterodimers_protein <- df_Xie2025_heterodimers_protein %>%
  rename("TF2" = "symbol")

Xie2025_protein_sequences <- protein_sequences




#How many proteins we get

nrow(df_Jolma2013_protein)+
nrow(df_Nitta2015)+nrow(df_motif_info %>% filter(study=="Morgunova2015"))+
nrow(df_Jolma2015_heterodimers_protein)+
nrow(df_Jolma2015_monomers_protein)+
nrow(df_Xie2025_monomers_protein)+
nrow(df_Xie2025_heterodimers_protein)+
nrow(df_Yin2017_protein) # 3636

tested=df_motif_info %>% filter(!(experiment %in% "Methyl-HT-SELEX")) # 3636

names(df_Jolma2013_protein)
names(df_Nitta2015)
names(df_Yin2017_protein)
df_Jolma2015_monomers_protein=df_Jolma2015_monomers_protein %>% rename(TF1=symbol)
df_Xie2025_monomers_protein =df_Xie2025_monomers_protein %>% rename(TF1=symbol)


df_Jolma2013_protein$symbol_clone=NULL
df_Nitta2015$symbol_clone=NULL
df_Yin2017_protein$symbol_clone=NULL
df_Jolma2015_monomers_protein$symbol_clone=NULL
df_Xie2025_monomers_protein$symbol_clone=NULL


 all_data=bind_rows(df_Jolma2013_protein,
  df_Nitta2015,
  df_motif_info %>% filter(study=="Morgunova2015"),
  df_Jolma2015_monomers_protein,
  df_Xie2025_monomers_protein,
  df_Yin2017_protein,
  df_Jolma2015_heterodimers_protein,
  df_Xie2025_heterodimers_protein
  )

all_data$`TF2 Classificaton in Vaquerizas et al., 2009`=NA
all_data$`TF2 production method`=NA
all_data$TF2_alternative_Ensembl_ID=NA
all_data$`TF2 Vaquerizas_comment`=NA


names(all_data)

all_data=all_data[,c("ID","old_ID", "motif_ID", "symbol","symbol_both", "TF1", "TF2", "Human_Ensemble_ID", "clone",
  "family",                                       "Lambert2018_families",
  "organism",                                     "study",                                        "experiment",
  "ligand",                                       "batch",                                        "seed",
  "multinomial",                                  "cycle",                                        "representative",
  "short",                                        "type",                                         "comment",
  "filename",                                     "IC",                                           "length",
  "consensus",                                    "kld_between_revcomp",                          "Methyl.SELEX.Motif.Category",
  "match_number",                                 "threshold",                                    "phyloP_threshold",
  "TF1 HNGC",                                     "TF1 Ensembl_ID",
  "TF1_protein_sequence",                         "TF1 DOMAIN",                                   "TF1 6-mer",
  "TF1 source organism",                          "TF1 Clone source",                             "TF1 Construct type",
  "TF1 Protein amount",                           "TF1 HT-SELEX validation",                      "TF1 Activity in CAP-SELEX",
  "TF1 Included domain Ids",                      "TF1 Included domain description",              "TF1 Excluded domain Ids",
  "TF1 Excluded domain description",              "TF1 clone type",
  "TF1 Classificaton in Vaquerizas et al., 2009", "TF1 production method",
  "TF1_alternative_Ensembl_ID",                   "TF1 Vaquerizas_comment",
  "TF2 HNGC",                                     "TF2 Ensembl_ID",
  "TF2_protein_sequence",                         "TF2 DOMAIN",                                   "TF2 6-mer",
  "TF2 source organism",                          "TF2 Clone source",                             "TF2 Construct type",
  "TF2 Protein amount",                           "TF2 HT-SELEX validation",                      "TF2 Activity in CAP-SELEX",
  "TF2 Included domain Ids",                      "TF2 Included domain description",              "TF2 Excluded domain Ids",
  "TF2 Excluded domain description",              "TF2 clone type",
  "TF2 Classificaton in Vaquerizas et al., 2009", "TF2 production method",
  "TF2_alternative_Ensembl_ID",                   "TF2 Vaquerizas_comment"
  )  ]

all_data$symbol %>% table(useNA = "always")
all_data$TF1 %>% table(useNA = "always")
all_data$TF2 %>% table(useNA = "always")

all_data$TF1[is.na(all_data$TF1)]=all_data$symbol[is.na(all_data$TF1)]

save(all_data, df_Jolma2013_protein,
  Jolma2013_protein_sequences,
  Jolma2013_protein_sequences_unique,
  df_Yin2017_protein,
  df_Yin2017_Methyl_protein,
  Yin2017_protein_sequences,
  Yin2017_protein_sequences_unique,
  df_Jolma2015_monomers_protein,
  df_Jolma2015_heterodimers_protein,
  Jolma2015_protein_sequences,
  df_Xie2025_monomers_protein,
  df_Xie2025_heterodimers_protein,
  Xie2025_protein_sequences,
  file = "~/projects/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData"
)

#For how many proteins we found protein sequences

#Monomers
all_data %>% filter(experiment!="CAP-SELEX") %>% nrow() # 1738
all_data %>% filter(experiment!="CAP-SELEX") %>% filter(!is.na(TF1_protein_sequence)) %>% nrow() # 1726 of 1738

test=all_data %>% filter(experiment!="CAP-SELEX") %>% filter(is.na(TF1_protein_sequence))
test$study %>% table()
#Morgunova2015     Nitta2015       Yin2017
#            1            10             1
#Heterodimers

all_data %>% filter(experiment=="CAP-SELEX") %>% nrow() # 1898

all_data %>% filter(experiment=="CAP-SELEX") %>% filter(!is.na(TF1_protein_sequence) & !is.na(TF2_protein_sequence)) %>% nrow() # 1897

test=all_data %>% filter(experiment=="CAP-SELEX") %>% filter(is.na(TF1_protein_sequence) | is.na(TF2_protein_sequence))
test$study %>% table()
#Xie2025
#      1
