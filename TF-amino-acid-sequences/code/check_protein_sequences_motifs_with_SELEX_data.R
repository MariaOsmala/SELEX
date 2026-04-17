#region packages

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(httr2)
library(tibble)
library(biomaRt)
library(Biostrings)
library(httr)

#regionend


#region load data
rm(list = ls())
data_path <- "/Users/osmalama/projects/SELEX/motif-metadata/"

data_file <- "motifs_with_SELEX_signal_and_background_and_proteins_04032026.tsv"

metadata <- read_delim(paste0(data_path,  data_file),
  col_types = cols(.default = col_character()),
 delim = "\t") #3594

metadata <- metadata |> filter(!is.na(`protein sequence 1`)) #3225

metadata <- metadata |>
  mutate(
    TF1_protein_sequence = `protein sequence 1`,
    TF2_protein_sequence =
      case_when(
        experiment == "CAP-SELEX" &
        TF1_copies == 1 ~ `protein sequence 2`, # nolint
        experiment == "CAP-SELEX" &
        TF1_copies == 2 ~ `protein sequence 3`,# nolint
        experiment == "CAP-SELEX" &
        TF1_copies == 3 ~ `protein sequence 4`,# nolint
        TRUE ~ NA_character_
      )
  )


TF_Info_Path <- "~/projects/TFBS/code/TF_interaction_summary/Data/" # nolint

TF_info_Ensembl115=read_delim( # nolint
  paste0(TF_Info_Path, # nolint.
         "TF_info_Ensembl115.txt"),
  delim = "\t"
)

length(unique(c(metadata$TF1, metadata$TF2))) #683 + NA

table(
  unique(c(metadata$TF1, metadata$TF2)) %in% TF_info_Ensembl115$symbol
) # all 683 are in the TF_info_Ensembl115

#endregion

#region build table
# Expected columns: motif, hgnc_symbol1, hgnc_symbol2, 
# organism1, organism2, ensembl_gene_id1, ensembl_gene_id1, 
# dbd_aa1, dbd_aa2

#Remove problematic cases

metadata <- metadata[-which(metadata$TF1 == "CART1"),] #3224
metadata$TF2[which(metadata$TF2 == "CART1")] #None


metadata <- metadata |> left_join(TF_info_Ensembl115 |> dplyr::select(
  symbol,
  organism,
  hgnc_symbol,
  symbol_is_hgnc_symbol,
  external_gene_name,
  description,
  ensembl_gene_id1,
), by = c("TF1" = "symbol")) |>
  dplyr::rename(hgnc_symbol1 = hgnc_symbol,
    organism1 = organism,
    TF1_ensembl_gene_id = ensembl_gene_id1,
    TF1_hgnc_symbol = hgnc_symbol,
    TF1_symbol_is_hgnc_symbol = symbol_is_hgnc_symbol,
    TF1_external_gene_name = external_gene_name,
    TF1_description = description
  )


metadata <- metadata |> left_join(TF_info_Ensembl115 |> dplyr::select(
  symbol,
  organism,
  hgnc_symbol,
  symbol_is_hgnc_symbol,
  external_gene_name,
  description,
  ensembl_gene_id1,
), by = c("TF2" = "symbol")) |>
  dplyr::rename(hgnc_symbol2 = hgnc_symbol,
    organism2 = organism,
    TF2_ensembl_gene_id = ensembl_gene_id1,
    TF2_hgnc_symbol = hgnc_symbol,
    TF2_symbol_is_hgnc_symbol = symbol_is_hgnc_symbol,
    TF2_external_gene_name = external_gene_name,
    TF2_description = description
  )

#regionend

#region functions
clean_aa <- function(x) {
  x <- toupper(ifelse(is.na(x), "", x))
  x <- gsub("\\s+", "", x)
  gsub("[^ACDEFGHIKLMNPQRSTVWYBXZJUO]", "", x)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#regionend

#region Ensembl
# 1) connect to Ensembl BioMart

httr_config <- switch(Sys.info()["sysname"],
                      "Linux" = config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
                      config())
set_config(config = httr_config)

options(timeout = 30000)
biomartCacheClear()
# Set up the Ensembl biomart for human genes
listEnsembl()$version[1]
# "Ensembl Genes 115" mirror = "useast"
mart <- biomaRt::useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 115
)


# 2) Get all peptide IDs per ENSG (i.e., all protein translations / isoforms)

# TF1

metadata_human <- metadata |>
  filter(organism1 != "Mus musculus")

ensg <- c(metadata_human$TF1_ensembl_gene_id,
          metadata_human$TF2_ensembl_gene_id) |>
  unique(na.rm = TRUE)



tx <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_peptide_id"
  ),
  filters    = "ensembl_gene_id",
  values     = ensg,
  mart       = mart
) |>
  filter(!is.na(ensembl_peptide_id), ensembl_peptide_id != "")

# 3) Fetch peptide sequences for those peptide IDs
#    biomaRt getSequence supports seqType="peptide" for protein sequences
pep_seqs <- getSequence(
  id      = unique(tx$ensembl_peptide_id),
  type    = "ensembl_peptide_id",
  seqType = "peptide",
  mart    = mart,
  useCache = TRUE
)

tx_seqs <- tx |>
  left_join(pep_seqs, by = "ensembl_peptide_id") |>
  filter(!is.na(peptide), peptide != "")

# Build a lookup: ENSG -> tibble(pep_id, peptide_seq)
by_gene <- split(tx_seqs, tx_seqs$ensembl_gene_id)

#region matching functions

# Exact containment (fast)
exact_hit <- function(dbd, peptide_vec) {
  if (dbd == "" || length(peptide_vec) == 0) return(NA_integer_)
  hits <- vmatchPattern(AAString(dbd), AAStringSet(peptide_vec), fixed = TRUE)
  idx <- which(length(hits) > 0)
  if (length(idx) == 0) NA_integer_ else idx[[1]]
}

# Strong alignment (slower, but robust to small differences)
# global-local = align whole DBD to a subsequence of the protein

best_alignment <- function(dbd, peptide_vec) {
  if (dbd == "" || length(peptide_vec) == 0) return(
    list(best_i = NA_integer_, pid2 = NA_real_, score = NA_real_)
  )
  alns <- lapply(peptide_vec, function(p) {
    pwalign::pairwiseAlignment(
      AAString(dbd),
      AAString(p),
      type = "global-local"
    )
  }
  )

  pid2  <- vapply(alns, function(a) pid(a, type = "PID2"), numeric(1))
  scr   <- vapply(alns, score, numeric(1))
  best  <- which.max(pid2)
  list(best_i = best, pid2 = pid2[[best]], score = scr[[best]])
}

#endregion

#region check sequences
# ---- row-wise check against claimed ENSG ----
metadata_checked <- metadata |>
  mutate(TF1_dbd_clean = clean_aa(TF1_protein_sequence)) |>
  rowwise() |>
  mutate(
    TF1_n_translations = nrow(by_gene[[TF1_ensembl_gene_id]] %||% tibble()),
    TF1_exact_translation_hit = {
      g <- by_gene[[TF1_ensembl_gene_id]]
      if (is.null(g)) {
        NA_character_
      } else {
        i <- exact_hit(TF1_dbd_clean, g$peptide)
        if (is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
      }
    },
    # Only align if exact containment fails:
    TF1_best_align = list({
      g <- by_gene[[TF1_ensembl_gene_id]]
      if (is.null(g) || !is.na(TF1_exact_translation_hit)) {
        list(best_i = NA_integer_, pid2 = NA_real_, score = NA_real_)
      } else {
        best_alignment(TF1_dbd_clean, g$peptide)
      }
    }),
    TF1_best_align_pid2 = TF1_best_align$pid2,
    TF1_best_align_score = TF1_best_align$score,
    TF1_best_align_translation_hit = {
      g <- by_gene[[TF1_ensembl_gene_id]]
      i <- TF1_best_align$best_i
      if (is.null(g) || is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
    },
    # choose a threshold you like; for many DBDs you might want >= 90–95% PID2
    TF1_aligns_strongly = is.na(TF1_exact_translation_hit) &&
      !is.na(TF1_best_align_pid2) && TF1_best_align_pid2 >= 90,
    TF1_dbd_matches_claimed_gene = !is.na(TF1_exact_translation_hit) ||
      TF1_aligns_strongly
  ) |>
  ungroup()

metadata_checked |>
  pull(TF1_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#FALSE  TRUE  <NA>
#  117  3107     0

metadata_checked |>
  filter(!TF1_dbd_matches_claimed_gene) |>
  dplyr::select(organism1) |>
  table(useNA = "ifany")
#Homo sapiens Mus musculus
#           1          116
#Failed for one Homo sapiens case, MYCL2, blasting gives the correct protein
metadata_checked |>
  filter(!TF1_dbd_matches_claimed_gene & organism1 == "Homo sapiens") |>
  dplyr::select(TF1, TF1_ensembl_gene_id, TF1_protein_sequence)


final <- metadata_checked |>
  mutate(TF2_dbd_clean = clean_aa(TF2_protein_sequence)) |>
  rowwise() |>
  mutate(
    TF2_n_translations = nrow(by_gene[[TF2_ensembl_gene_id]] %||% tibble()),
    TF2_exact_translation_hit = {
      g <- by_gene[[TF2_ensembl_gene_id]]

      if (is.null(g)) {
        NA_character_
      } else {
        i <- exact_hit(TF2_dbd_clean, g$peptide)
        if (is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
      }
    },
    # Only align if exact containment fails:
    TF2_best_align = list({
      g <- by_gene[[TF2_ensembl_gene_id]]
      if (is.null(g) || !is.na(TF2_exact_translation_hit)) {
        list(best_i = NA_integer_, pid2 = NA_real_, score = NA_real_)
      } else {
        best_alignment(TF2_dbd_clean, g$peptide)
      }
    }),
    TF2_best_align_pid2 = TF2_best_align$pid2,
    TF2_best_align_score = TF2_best_align$score,
    TF2_best_align_translation_hit = {
      g <- by_gene[[TF2_ensembl_gene_id]]
      i <- TF2_best_align$best_i
      if (is.null(g) || is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
    },
    # choose a threshold you like; for many DBDs you might want >= 90–95% PID2
    TF2_aligns_strongly = is.na(TF2_exact_translation_hit) &&
      !is.na(TF2_best_align_pid2) && TF2_best_align_pid2 >= 90,
    TF2_dbd_matches_claimed_gene = !is.na(TF2_exact_translation_hit) ||
      TF2_aligns_strongly
  ) |>
  ungroup()


final |>
  filter(!is.na(TF2)) |> nrow() #1640

final |>
  pull(TF2_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#FALSE  TRUE  <NA>
# 1584  1640     0

metadata <- final

#endregion

#region Mouse

biomartCacheClear()

options(timeout = 30000)
# Mouse
mart <- biomaRt::useEnsembl(biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  version = 115
) #mirror = "useast"


metadata_mouse <- metadata |>
  filter(organism1 == "Mus musculus")

ensg <- c(metadata_mouse$TF1_ensembl_gene_id) |> unique(na.rm = TRUE)

tx <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_peptide_id"
  ),
  filters    = "ensembl_gene_id",
  values     = ensg,
  mart       = mart
) |>
  filter(!is.na(ensembl_peptide_id), ensembl_peptide_id != "")

# 3) Fetch peptide sequences for those peptide IDs
#    biomaRt getSequence supports seqType="peptide" for protein sequences
pep_seqs <- getSequence(
  id      = unique(tx$ensembl_peptide_id),
  type    = "ensembl_peptide_id",
  seqType = "peptide",
  mart    = mart,
  useCache = TRUE
)

tx_seqs <- tx |>
  left_join(pep_seqs, by = "ensembl_peptide_id") |>
  filter(!is.na(peptide), peptide != "")

# Build a lookup: ENSG -> tibble(pep_id, peptide_seq)
by_gene <- split(tx_seqs, tx_seqs$ensembl_gene_id)

#endregion

#region check sequences

#does this overwrite the human checks? Yes, but we can re-run the human checks if needed. 
#We could also split the metadata into human and mouse parts and then combine after checking, but this is simpler for now.

# do not overwrite if the values exist

metadata_checked <- metadata_mouse |>
  mutate(TF1_dbd_clean = clean_aa(TF1_protein_sequence)) |>
  rowwise() |>
  mutate(
    TF1_n_translations = nrow(by_gene[[TF1_ensembl_gene_id]] %||% tibble()),
    TF1_exact_translation_hit = {
      g <- by_gene[[TF1_ensembl_gene_id]]
      if (is.null(g)) {
        NA_character_
      } else {
        i <- exact_hit(TF1_dbd_clean, g$peptide)
        if (is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
      }
    },
    # Only align if exact containment fails:
    TF1_best_align = list({
      g <- by_gene[[TF1_ensembl_gene_id]]
      if (is.null(g) || !is.na(TF1_exact_translation_hit)) {
        list(best_i = NA_integer_, pid2 = NA_real_, score = NA_real_)
      } else {
        best_alignment(TF1_dbd_clean, g$peptide)
      }
    }),
    TF1_best_align_pid2 = TF1_best_align$pid2,
    TF1_best_align_score = TF1_best_align$score,
    TF1_best_align_translation_hit = {
      g <- by_gene[[TF1_ensembl_gene_id]]
      i <- TF1_best_align$best_i
      if (is.null(g) || is.na(i)) NA_character_ else g$ensembl_peptide_id[[i]]
    },
    # choose a threshold you like; for many DBDs you might want >= 90–95% PID2
    TF1_aligns_strongly = is.na(TF1_exact_translation_hit) &&
      !is.na(TF1_best_align_pid2) && TF1_best_align_pid2 >= 90,
    TF1_dbd_matches_claimed_gene = !is.na(TF1_exact_translation_hit) ||
      TF1_aligns_strongly
  ) |>
  ungroup()

metadata_checked |>
  pull(TF1_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#TRUE  <NA>
#  116     0

metadata_checked |>
  dplyr::select(organism1) |>
  table(useNA = "ifany")
#116 Mus musculus motifs

rbind(final |> filter(organism1 != "Mus musculus"), metadata_checked) |>
  pull(TF1_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#FALSE  TRUE  <NA>
#    1(MYCL2)  3223     0

metadata <- rbind(final |> filter(organism1 != "Mus musculus"),
                  metadata_checked) #should be 3224, CART was removed

data_file <-
  "motifs_with_SELEX_signal_and_background_and_proteins_04032026_checked.tsv"

write_delim(metadata,
  paste0(data_path,
    data_file
  ),
  delim = "\t"
)

#endregion