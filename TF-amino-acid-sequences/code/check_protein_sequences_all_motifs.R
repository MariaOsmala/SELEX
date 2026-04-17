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

load(file = "~/projects/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData")
#all_data
metadata <- all_data # 3636
rm(list=setdiff(ls(), "metadata"))

(metadata |> filter(experiment!="CAP-SELEX") |> dplyr::select(TF1_protein_sequence) %>% is.na()) |> table(useNA = "always")
#FALSE  TRUE  <NA>
# 1726    12     0
(metadata |> filter(experiment=="CAP-SELEX") |> dplyr::select(TF2_protein_sequence) %>% is.na()) |> table(useNA = "always")
#FALSE  TRUE  <NA>
# 1898     1     0
#TF_protein_sequence present for all CAP-SELEX cases

TF_Info_Path <- "Data/" # nolint

TF_info_Ensembl115=read_delim( # nolint
  paste0(TF_Info_Path, # nolint.
         "TF_info_Ensembl115.tsv"),
  delim = "\t"
) # 728

TFs=unique(c(metadata$TF1, metadata$TF2))
TFs=TFs[-which(is.na(TFs))] # 728

table(
  TFs %in% TF_info_Ensembl115$symbol
)

TFs[!TFs %in% TF_info_Ensembl115$symbol]
#which TFs are not in the TF_info_Ensembl115?
#"Egr1_E410D_FARSDERtoFARSDDR"
TF_info_Ensembl115$symbol[which(TF_info_Ensembl115$symbol== "Egr1")]="Egr1_E410D_FARSDERtoFARSDDR"
table(
  TFs %in% TF_info_Ensembl115$symbol
) # 728

#Does TF_info_Ensembl115 contain all TFs with motifs, yes!
table(
  TF_info_Ensembl115$symbol %in% TFs
)


#endregion

#region build table
# Expected columns: motif, hgnc_symbol1, hgnc_symbol2,
# organism1, organism2, ensembl_gene_id1, ensembl_gene_id1,
# dbd_aa1, dbd_aa2
metadata$organism=NULL

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


biomartCacheClear()
options(timeout = 6000)
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

ensg=ensg[-which(is.na(ensg))]

#remove e coli TF
ensg=ensg[-which(ensg=="b3438")] # 637

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
  mart    = mart
)

tx_seqs <- tx |>
  left_join(pep_seqs, by = "ensembl_peptide_id") |>
  filter(!is.na(peptide), peptide != "")

# Build a lookup: ENSG -> tibble(pep_id, peptide_seq)
by_gene <- split(tx_seqs, tx_seqs$ensembl_gene_id)

#region matching functions

# Exact containment (fast)
exact_hit <- function(dbd, peptide_vec) {
  # exact_hit(TF1_dbd_clean, g$peptide)
  # dbd=TF1_dbd_clean
  # peptide_vec=g$peptide
  if (dbd == "" || length(peptide_vec) == 0) return(NA_integer_)
  hits <- vmatchPattern(AAString(dbd), AAStringSet(peptide_vec), fixed = TRUE)
  idx <- which(sapply(hits, length) > 0)
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

  pid2  <- vapply(alns, function(a) pwalign::pid(a, type = "PID2"), numeric(1))
  scr   <- vapply(alns, pwalign::score, numeric(1))
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
# 163  3473     0

metadata_checked |>
  filter(!TF1_dbd_matches_claimed_gene) |>
  dplyr::select(organism1) |>
  table(useNA = "ifany")
#escherichia_coli_str_k_12_substr_mg1655_gca_000005845   Homo sapiens  Mus musculus
#                                                    2   28            133

#Failed for one Homo sapiens case, MYCL2, blasting gives the correct protein
#15 problematic cases, alignment pid smaller than 70
problems=metadata_checked |>
  filter(!TF1_dbd_matches_claimed_gene & organism1 == "Homo sapiens" & !is.na(TF1_protein_sequence)) |>
  dplyr::select(ID,study,TF1, TF1_ensembl_gene_id, TF1_protein_sequence, TF1_best_align_pid2)

#Check which proteins these match
problems=metadata_checked |>
  filter(!TF1_dbd_matches_claimed_gene & organism1 == "Homo sapiens" & !is.na(TF1_protein_sequence)) |>
  dplyr::select(ID,TF1_protein_sequence)

# write one multi-FASTA file
fasta_lines <- as.vector(rbind(paste0(">", problems$ID), problems$TF1_protein_sequence))
writeLines(fasta_lines, "data/missing_human_proteins.fa")

#do blastp agains uniprot

blastp=read_delim("Data/Missing_human_proteins_blast.csv",delim = ",", col_names=FALSE)
head(blastp)
names(blastp)=c("ID", "hit", "pid","alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "% positives")
blastp <- split(blastp, blastp$ID)
blastp<- do.call(rbind, lapply(blastp, function(x) x[which(x$pid == max(x$pid)), ]))
blastp$uniprot_gn_id=do.call(rbind, strsplit(blastp$hit,"\\."))[,1]
res <- getBM(
  attributes = c("uniprot_gn_symbol" ,"uniprot_gn_id", "ensembl_gene_id", "hgnc_symbol", "description"),
  filters    = "uniprot_gn_id",
  values     = blastp$uniprot_gn_id,
  mart       = mart
)

blastp <- blastp |> left_join(res, by = "uniprot_gn_id")
blastp <- blastp%>% relocate(c(uniprot_gn_symbol, hgnc_symbol, description), .before = hit)

metadata_checked=metadata_checked |>
  mutate(TF1_blastp_hit_uniprot_gn_id = blastp$uniprot_gn_id[match(ID, blastp$ID)],
          TF1_blastp_hit_ensembl_ID = blastp$ensembl_gene_id[match(ID, blastp$ID)],
          TF1_blastp_hit_hgnc_symbol = blastp$hgnc_symbol[match(ID, blastp$ID)],
          TF1_blastp_hit_description = blastp$description[match(ID, blastp$ID)],
          TF1_blastp_hit_pid = blastp$pid[match(ID, blastp$ID)]
         )

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
  filter(!is.na(TF2)) |> nrow() # 1898

final |>
  pull(TF2_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#FALSE  TRUE  <NA>
# 1738  1898     0

metadata <- final

#endregion

#region Mouse

biomartCacheClear()

options(timeout = 3000)
# Mouse
mart <- biomaRt::useEnsembl(biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  version = 115
) #mirror = "useast"

metadata_mouse <- metadata |>
  filter(organism1 == "Mus musculus")

ensg <- c(metadata_mouse$TF1_ensembl_gene_id) |> unique(na.rm = TRUE)
is.na(ensg) %>% table(useNA = "always")

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
  mart    = mart
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
#  133     0

metadata_checked |>
  dplyr::select(organism1) |>
  table(useNA = "ifany")
# 133 Mus musculus motifs

rbind(final |> filter(organism1 != "Mus musculus"), metadata_checked) |>
  pull(TF1_dbd_matches_claimed_gene) |>
  table(useNA = "always")
#FALSE  TRUE  <NA>
#    30 3606     0

metadata <- rbind(final |> filter(organism1 != "Mus musculus"),
                  metadata_checked) #should be 3636ed

#endregion

#region E colil
#E coli TF, which is not in Ensembl, but we can check the sequence manually:

metadata_ecoli <- metadata |>
  filter(organism1 == "escherichia_coli_str_k_12_substr_mg1655_gca_000005845")

ensg <- c(metadata_ecoli$TF1_ensembl_gene_id) |> unique(na.rm = TRUE)



library(jsonlite)
library(dplyr)

# 2) Get all linked xrefs for that gene, including lower levels
url2 <- sprintf(
  "https://rest.ensembl.org/xrefs/id/%s?all_levels=1;content-type=application/json",
  ensg
)

xrefs <- fromJSON(url2) |> as_tibble()
xrefs

#Next, find the linked protein/translation identifier. Depending on the record, this may already be in the symbol
# lookup output or may need to be taken from the gene lookup:

# 3) Also do a lookup on the gene ID to inspect its structure
url3 <- sprintf(
  "https://rest.ensembl.org/lookup/id/%s?expand=1;content-type=application/json",
  ensg
)

gene_lookup <- fromJSON(url3)
gene_lookup

#For bacterial genes, there is usually one transcript and one translation.
#Extract the translation ID and fetch the protein sequence:

# 4) Pull translation/protein Ensembl ID from expanded lookup
translation_id <- gene_lookup$Transcript$Translation$id
translation_id

# 5) Get the protein sequence
url4 <- sprintf(
  "https://rest.ensembl.org/sequence/id/%s?type=protein;content-type=text/plain",
  translation_id
)

ensembl_protein_seq <- readLines(url4, warn = FALSE) |> paste(collapse = "")
ensembl_protein_seq

aa_query=metadata |>
  filter(organism1 == "escherichia_coli_str_k_12_substr_mg1655_gca_000005845") |>
  pull(TF1_protein_sequence) %>% unique()
# MKKKRPVLQDVADRVGVTKMTVSRFLRNPEQVSVALRGKIAAALDELGYIPNRAPDILSNATSRAIGVLLPSLTN
# 6) Compare
identical(aa_query, ensembl_protein_seq) #FALSE

nchar(aa_query)
nchar(ensembl_protein_seq)

# Simple prefix / containment checks
startsWith(ensembl_protein_seq, aa_query) #TRUE
grepl(aa_query, ensembl_protein_seq, fixed = TRUE) #TRUE

metadata <- metadata %>%
  mutate("TF1_dbd_matches_claimed_gene" = if_else(organism1 == "escherichia_coli_str_k_12_substr_mg1655_gca_000005845",
   TRUE,
   TF1_dbd_matches_claimed_gene)
  )

metadata %>% dplyr::select(experiment) %>% table(useNA = "always")
metadata %>% filter(experiment == "HT-SELEX") %>% dplyr::select(TF1_dbd_matches_claimed_gene) %>% table(useNA = "always")
#FALSE  TRUE  <NA>
#   26  1712     0
missing_HTSELEX=metadata %>% filter(experiment == "HT-SELEX" & TF1_dbd_matches_claimed_gene == FALSE) %>% dplyr::select(TF1,TF1_protein_sequence,TF1_dbd_matches_claimed_gene)
#MYCL1 or protein missing

metadata %>% filter(experiment == "CAP-SELEX" & TF1_dbd_matches_claimed_gene == FALSE) %>% dplyr::select(TF1,TF1_protein_sequence,TF1_dbd_matches_claimed_gene)
#ALX3 and TCF15 missing
metadata %>% filter(experiment == "CAP-SELEX" & TF2_dbd_matches_claimed_gene == FALSE) %>% dplyr::select(TF2,TF2_protein_sequence,TF2_dbd_matches_claimed_gene)
#Empty

metadata %>% dplyr::select(TF1_dbd_matches_claimed_gene) %>% table(useNA = "always")
#FALSE  TRUE  <NA>
#   28  3608     0

data_file <-
  "motifs_with_proteins.tsv"

write_delim(metadata,
  paste0("motif-metadata/",
    data_file
  ),
  delim = "\t"
)

#endregion
