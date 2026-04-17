library(httr)
library("readr")
# Load the biomaRt library
library(biomaRt)
library("dplyr")
library("stringr")
library("tidyverse")


rm(list=ls())

TFs_with_motifs=readRDS(file="RData/TFs_with_motifs.RDS") # 728

########## GRCh38 from Ensemble ##################

httr_config <- switch(Sys.info()["sysname"],
                      "Linux" = config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
                      config())
set_config(config = httr_config)

biomartCacheClear()
options(timeout = 600)  # 10 minutes timeout
# Set up the Ensembl biomart for human genes
#listEnsembl()$version[1]
#"Ensembl Genes 115"
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version=115) #mirror = "useast"

message("Connection to Ensembl established")

# Query all human genes, 795865



#length(TFs_with_motifs) # 728

has_lower <- which(str_detect(TFs_with_motifs, "\\p{Ll}") ) #83
mouse_TFs_with_motifs=data.frame(symbol=TFs_with_motifs[has_lower])
human_TFs_with_motifs=data.frame(symbol=TFs_with_motifs[-has_lower]) # 645
ecoli_TFs_with_motifs=data.frame(symbol= mouse_TFs_with_motifs[nrow(mouse_TFs_with_motifs),]) # 83
mouse_TFs_with_motifs=data.frame(symbol= mouse_TFs_with_motifs[-nrow(mouse_TFs_with_motifs),]) # 0

filters <- list("hgnc_symbol"=human_TFs_with_motifs$symbol)
attributes <- c("hgnc_symbol","ensembl_gene_id", "external_synonym","external_gene_name", "entrezgene_accession", "description")
result <- getBM(attributes, filters, mart = mart)

df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(hgnc_symbol, {{col}})

  if (sort_vals) out <- out %>% arrange(hgnc_symbol, {{col}})

  out %>%
    group_by(hgnc_symbol) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(hgnc_symbol) %>%
  summarise(
    external_gene_name   = dplyr::first(na.omit(external_gene_name)),
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )



# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
entrez_wide <- wide_many(df, entrezgene_accession, max_n = 10)
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, entrez_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "hgnc_symbol")

final_wide

human_TFs_with_motifs %>% nrow() # 645

human_TFs_with_motifs=left_join(human_TFs_with_motifs, final_wide, by=c("symbol"="hgnc_symbol"), keep=TRUE)

#For some TFs, no info was found with the "hgnc_symbol" filter,

#try with "external_synonym"  or "external_gene_name" (did not work)

filters <- list("external_synonym" = human_TFs_with_motifs %>% dplyr::filter(is.na(hgnc_symbol)) %>% pull(symbol)) #15
attributes <- c("hgnc_symbol", "ensembl_gene_id","external_synonym","external_gene_name", "entrezgene_accession","description")
result=getBM(attributes, filters, mart = mart)


df <- result %>% distinct()   # de-dup exact rows, just in case


# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(hgnc_symbol, {{col}})

  if (sort_vals) out <- out %>% arrange(hgnc_symbol, {{col}})

  out %>%
    group_by(hgnc_symbol) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(hgnc_symbol) %>%
  summarise(
    external_gene_name   = dplyr::first(na.omit(external_gene_name)),
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
entrez_wide <- wide_many(df, entrezgene_accession, max_n = 10)
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, entrez_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "hgnc_symbol")

final_wide


human_TFs_with_motifs=left_join(human_TFs_with_motifs, final_wide, by=c("symbol"="external_synonym1"), keep=TRUE)

fill_ind=which(!is.na(human_TFs_with_motifs$hgnc_symbol.y ))

human_TFs_with_motifs[fill_ind, c("hgnc_symbol.x",
                 "external_gene_name.x",
                 "description.x",
                 "ensembl_gene_id1.x",
                 "entrezgene_accession1.x",
                 "external_synonym1.x"  )]= human_TFs_with_motifs[fill_ind,
                  c("hgnc_symbol.y",
                  "external_gene_name.y",
                  "description.y",
                  "ensembl_gene_id1.y",
                  "entrezgene_accession1.y",
                  "external_synonym1.y")]

#remove columns that end with .y

human_TFs_with_motifs=human_TFs_with_motifs %>% dplyr::select(-ends_with(".y"))

#remove .x from column names

names(human_TFs_with_motifs)=gsub("\\.x$", "", names(human_TFs_with_motifs))

#Now contains 646 rows instead of 645, CART1 is twice
#CART1 is likely ALX1, as CART1 clusters together with ALX1 motifs
human_TFs_with_motifs %>% dplyr::filter(symbol=="CART1")
# symbol hgnc_symbol external_gene_name                                                          description ensembl_gene_id1 ensembl_gene_id2 ensembl_gene_id3
# 1  CART1        ALX1               ALX1                    ALX homeobox 1 [Source:HGNC Symbol;Acc:HGNC:1494]  ENSG00000180318             <NA>             <NA>
# 2  CART1       TRAF4              TRAF4 TNF receptor associated factor 4 [Source:HGNC Symbol;Acc:HGNC:12034]  ENSG00000076604

# Remove one of the CART1 rows

human_TFs_with_motifs=human_TFs_with_motifs %>% dplyr::filter(!(symbol=="CART1" & hgnc_symbol=="TRAF4")) # 645

#still something missing

human_TFs_with_motifs %>% dplyr::filter(is.na(hgnc_symbol)) %>% pull(symbol) #"HINFP1", this is likely HINFP

filters <- list("hgnc_symbol"="HINFP")
attributes <- c("hgnc_symbol","ensembl_gene_id", "external_synonym","external_gene_name", "entrezgene_accession", "description")
result <- getBM(attributes, filters, mart = mart)

df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(hgnc_symbol, {{col}})

  if (sort_vals) out <- out %>% arrange(hgnc_symbol, {{col}})

  out %>%
    group_by(hgnc_symbol) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(hgnc_symbol) %>%
  summarise(
    external_gene_name   = dplyr::first(na.omit(external_gene_name)),
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
entrez_wide <- wide_many(df, entrezgene_accession, max_n = 10)
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, entrez_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "hgnc_symbol")

final_wide

missing_ind=which(human_TFs_with_motifs$symbol=="HINFP1")


human_TFs_with_motifs[missing_ind, c("hgnc_symbol",
                 "external_gene_name",
                 "description",
                 "ensembl_gene_id1",
                 "entrezgene_accession1",
                 "external_synonym1"  )]= final_wide[1,
                  c("hgnc_symbol",
                  "external_gene_name",
                  "description",
                  "ensembl_gene_id1",
                  "entrezgene_accession1",
                  "external_synonym1")]

table(is.na(human_TFs_with_motifs$ensembl_gene_id1)) #All false

# Uniprot IDs -------------------------------------------------------------


# "ensembl_gene_id_version", , "uniprot_gn_id", "uniprot_gn_symbol" , "uniprot_gn_id", "uniprot_gn_symbol"

# 4) Query by a list of gene symbols
result <- getBM(
  attributes = c(
    "ensembl_gene_id", "uniprot_gn_id", "uniprot_gn_symbol"),
  filters   = list("ensembl_gene_id"=human_TFs_with_motifs$ensembl_gene_id1),
  mart      = mart
)

df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(ensembl_gene_id, {{col}})

  if (sort_vals) out <- out %>% arrange(ensembl_gene_id, {{col}})

  out %>%
    group_by(ensembl_gene_id) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    uniprot_gn_symbol         = dplyr::first(na.omit(uniprot_gn_symbol)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns

uniprot_wide    <- wide_many(df, uniprot_gn_id, max_n = 20)

# Final one-row-per-gene table
final_wide1 <- list(base, uniprot_wide) %>%
  purrr::reduce(left_join, by = "ensembl_gene_id")

final_wide1

#For some TFs, there are different symbol names

(human_TFs_with_motifs$symbol!=human_TFs_with_motifs$hgnc_symbol) %>% table()
# FALSE  TRUE
# 630    15

human_TFs_with_motifs$symbol_is_hgnc_symbol=(human_TFs_with_motifs$symbol==human_TFs_with_motifs$hgnc_symbol)

human_TFs_with_motifs =human_TFs_with_motifs %>% relocate(symbol_is_hgnc_symbol, .after=hgnc_symbol)
human_TFs_with_motifs=left_join(human_TFs_with_motifs, final_wide1, by=c("ensembl_gene_id1"="ensembl_gene_id"))

missing=human_TFs_with_motifs %>% filter(uniprot_gn_symbol=="") #no uniprot data for two genes


# Same for mouse TFs ------------------------------------------------------

nrow(mouse_TFs_with_motifs) #82
httr_config <- switch(Sys.info()["sysname"],
                      "Linux" = config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
                      config())
set_config(config = httr_config)
biomartCacheClear()
# Set up the Ensembl biomart for human genes
#listEnsembl()$version[1]
#"Ensembl Genes 115"
options(timeout = 600)  # 10 minutes timeout
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version=115) #mirror = "useast"

message("Connection to Ensembl established")

# 2) See what you can ask for
#   listAttributes(mart)   # huge table of available columns
#   listFilters(mart)      # possible filters
#   searchAttributes(mart, "homolog")  # quick search
#   searchAttributes(mart, "transcript")
# "ensembl_gene_id_version", , "uniprot_gn_id", "uniprot_gn_symbol" , "uniprot_gn_id", "uniprot_gn_symbol"

# 4) Query by a list of gene symbols
result <- getBM(
  attributes = c(
    "ensembl_gene_id", "external_synonym", "external_gene_name", "description"),
  filters   = list("external_gene_name"=mouse_TFs_with_motifs$symbol),
  mart      = mart
)

df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(external_gene_name, {{col}})

  if (sort_vals) out <- out %>% arrange(external_gene_name, {{col}})

  out %>%
    group_by(external_gene_name) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(external_gene_name) %>%
  summarise(
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "external_gene_name")

final_wide #78

mouse_TFs_with_motifs %>% nrow() #82

mouse_TFs_with_motifs=left_join(mouse_TFs_with_motifs, final_wide, by=c("symbol"="external_gene_name"), keep=TRUE)

#For some TFs, no info was found with the "external_gene_name" filter,

#try with "external_synonym"

filters <- list("external_synonym" = mouse_TFs_with_motifs %>% dplyr::filter(is.na(external_gene_name)) %>% pull(symbol)) #5
attributes <- c("ensembl_gene_id","external_synonym","external_gene_name", "description")
result=getBM(attributes, filters, mart = mart)

df <- result %>% distinct()   # de-dup exact rows, just in case


# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(external_gene_name, {{col}})

  if (sort_vals) out <- out %>% arrange(external_gene_name, {{col}})

  out %>%
    group_by(external_gene_name) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(external_gene_name) %>%
  summarise(
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "external_gene_name")

final_wide



mouse_TFs_with_motifs=left_join(mouse_TFs_with_motifs, final_wide, by=c("symbol"="external_synonym1"), keep=TRUE)

fill_ind=which(!is.na(mouse_TFs_with_motifs$external_gene_name.y ))

mouse_TFs_with_motifs[fill_ind, c(
                 "external_gene_name.x",
                 "description.x",
                 "ensembl_gene_id1.x",
                 "external_synonym1.x"  )]= mouse_TFs_with_motifs[fill_ind,
                  c(
                  "external_gene_name.y",
                  "description.y",
                  "ensembl_gene_id1.y",
                  "external_synonym1.y")]

#remove columns that end with .y

mouse_TFs_with_motifs=mouse_TFs_with_motifs %>% dplyr::select(-ends_with(".y"))

#remove .x from column names

names(mouse_TFs_with_motifs)=gsub("\\.x$", "", names(mouse_TFs_with_motifs))


#still something missing

mouse_TFs_with_motifs %>% dplyr::filter(is.na(external_gene_name)) %>% pull(symbol)
#"Tp53" "Tp73"
#Trp53 and Trp73 are mouse homologs of human TP53 and TP73, respectively


filters <- list("external_gene_name"=c("Trp53", "Trp73"))
attributes <- c("ensembl_gene_id", "external_synonym","external_gene_name", "description")
result <- getBM(attributes, filters, mart = mart)


df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(external_gene_name, {{col}})

  if (sort_vals) out <- out %>% arrange(external_gene_name, {{col}})

  out %>%
    group_by(external_gene_name) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(external_gene_name) %>%
  summarise(
    description          = dplyr::first(na.omit(description)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns
ens_wide    <- wide_many(df, ensembl_gene_id,      max_n = 10)  # change max_n if desired
syn_wide    <- wide_many(df, external_synonym,     max_n = 20)

# Final one-row-per-gene table
final_wide <- list(base, ens_wide, syn_wide) %>%
  purrr::reduce(left_join, by = "external_gene_name")

final_wide

missing_ind=match(gsub("r", "", final_wide$external_gene_name),mouse_TFs_with_motifs$symbol)



mouse_TFs_with_motifs[missing_ind, c(
                 "external_gene_name",
                 "description",
                 "ensembl_gene_id1",
                 "external_synonym1","external_synonym2","external_synonym3"  )]= final_wide[1:2,
                  c(
                  "external_gene_name",
                  "description",
                  "ensembl_gene_id1",
                  "external_synonym1",
                  "external_synonym2","external_synonym3")]

table(is.na(mouse_TFs_with_motifs$ensembl_gene_id1))
#Uniprot IDs of the TFs

# Uniprot IDs -------------------------------------------------------------


# "ensembl_gene_id_version", , "uniprot_gn_id", "uniprot_gn_symbol" , "uniprot_gn_id", "uniprot_gn_symbol"

# 4) Query by a list of gene symbols
result <- getBM(
  attributes = c(
    "ensembl_gene_id", "uniprot_gn_id", "uniprot_gn_symbol"),
  filters   = list("ensembl_gene_id"=mouse_TFs_with_motifs$ensembl_gene_id1),
  mart      = mart
)

df <- result %>% distinct()   # de-dup exact rows, just in case

# Helper: make a wide block for one multi-valued column
wide_many <- function(data, col, prefix = deparse(substitute(col)), max_n = Inf, sort_vals = FALSE) {
  out <- data %>%
    dplyr::filter(!is.na({{col}}) & {{col}} != "") %>%
    distinct(ensembl_gene_id, {{col}})

  if (sort_vals) out <- out %>% arrange(ensembl_gene_id, {{col}})

  out %>%
    group_by(ensembl_gene_id) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    dplyr::filter(idx <= max_n) %>%                                    # cap number of columns if you want
    pivot_wider(names_from = idx, values_from = {{col}},
                names_glue = paste0(prefix, "{idx}"))
}

# "Single-valued" info (pick a representative non-NA value per gene)
base <- df %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    uniprot_gn_symbol         = dplyr::first(na.omit(uniprot_gn_symbol)),
    .groups = "drop"
  )

# Wide blocks for the multi-valued columns

uniprot_wide    <- wide_many(df, uniprot_gn_id, max_n = 20)

# Final one-row-per-gene table
final_wide1 <- list(base, uniprot_wide) %>%
  purrr::reduce(left_join, by = "ensembl_gene_id")

final_wide1

mouse_TFs_with_motifs=left_join(mouse_TFs_with_motifs, final_wide1, by=c("ensembl_gene_id1"="ensembl_gene_id"))


names(mouse_TFs_with_motifs)
names(human_TFs_with_motifs)

human_TFs_with_motifs$organism="Homo sapiens"
mouse_TFs_with_motifs$organism="Mus musculus"

# E coli TF

library(jsonlite)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# 1) symbol -> Ensembl objects
species = "escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
url1 <- sprintf(
  "https://rest.ensembl.org/xrefs/symbol/%s/%s?content-type=application/json",
  species, utils::URLencode(ecoli_TFs_with_motifs$symbol, reserved = TRUE)
)

hits <- jsonlite::fromJSON(url1, simplifyDataFrame = TRUE)
hits <- tibble::as_tibble(hits)

# Prefer a gene-level hit
gene_hits <- hits %>% filter(type == "gene")
gene_id <- gene_hits$id[[1]]


# 2) Ensembl ID -> external refs
url2 <- sprintf(
  "https://rest.ensembl.org/xrefs/id/%s?all_levels=1;content-type=application/json",
  gene_id
)

xrefs <- jsonlite::fromJSON(url2, simplifyDataFrame = TRUE)
xrefs <- tibble::as_tibble(xrefs)

uniprot_id=xrefs %>% filter(db_display_name== "UniProtKB/Swiss-Prot")%>%pull( "display_id")

ecoli_TFs_with_motifs$ensembl_gene_id1=gene_id
ecoli_TFs_with_motifs$uniprot_gn_id1 = uniprot_id
ecoli_TFs_with_motifs$organism=species





names(human_TFs_with_motifs)
names(mouse_TFs_with_motifs)

combined_TFs <- bind_rows(
  human_TFs_with_motifs,
  mouse_TFs_with_motifs,
  ecoli_TFs_with_motifs
)

# relocate organism column

combined_TFs=combined_TFs %>% relocate(organism, .after=symbol)

write.table(combined_TFs, file="Data/TF_info_Ensembl115.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(combined_TFs, file="RData/TF_info_Ensembl115.RDS")


combined_TFs=readRDS(file="RData/TF_info_Ensembl115.RDS")

# Add also Lambert2018 info and TFClass info


tfclass=readRDS(file = "~/projects/TFBS/code/motif-networks/TFClass/tfclass.RDS")


#Try to add also Lambert2018 data here

Lambert2018 <- read_excel("~/projects/TFBS/code/motif-networks/TFClass/Lambert2018_supplementary_tables_S1_S4.xlsx",
                          sheet = "Table S1. Related to Figure 1B",
                          col_types="text",
                          skip = 1)

Lambert2018$Name_upper <- toupper(Lambert2018$Name)
Lambert2018$Name_upper_without_dash <- gsub("-", "",Lambert2018$Name_upper)

Lambert2018$Name_upper_without_dash %>% unique %>% length() #2765

#Try to match cytoscape network data to Lambert2018 data and TFClass data

monomers <- read_csv("~/projects/cytoscape/monomer_network/node_table_monomer_network_logos.csv")

monomers %>% select(node_type) %>% table()
#logo  motif    TF
#2035  2035   718

monomers %>% filter(node_type=="TF") %>% select(name) %>% head()
monomers %>% filter(node_type=="motif") %>% select(name) %>% head()

monomer_motifs=monomers %>% filter(node_type=="motif")

monomer_TFs=monomers %>% filter(node_type=="TF")

#Exctract the family info from monomer_motifs

match_inds=lapply(monomer_TFs$name, function(x) which(!is.na(match(monomer_motifs$symbol, x))))

orig_families=lapply(match_inds, function(x) unique(monomer_motifs$family[x]))

#Separate by ;

combined <- sapply(orig_families, function(x) {
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
})

monomer_TFs$family <- combined

#Remove NA; and NA

monomer_TFs$family <- gsub("NA; ", "", monomer_TFs$family)
monomer_TFs$family <- gsub("; NA", "", monomer_TFs$family)
monomer_TFs$family <- gsub("NA", "", monomer_TFs$family)

monomer_TFs$Name= toupper(monomer_TFs$name)

monomer_TFs$Name_without_dash= gsub("-", "", monomer_TFs$Name)

monomer_TFs$Name %in% Lambert2018$Name %>% table()

#FALSE  TRUE
#17   701

test=monomer_TFs %>% filter(!(Name %in% Lambert2018$Name))

monomer_TFs_Lambert2018 <- monomer_TFs %>%
  left_join(Lambert2018, by = "Name")

monomer_TFs_Lambert2018=rename(monomer_TFs_Lambert2018, Ensembl_ID = ID)
monomer_TFs_Lambert2018=rename(monomer_TFs_Lambert2018, `Is TF?` = "...4" )
monomer_TFs_Lambert2018=rename(monomer_TFs_Lambert2018, `Lambert2018 PDB` = PDB )

unique_tfclass_gene_upper_without_dash= unique(c(tfclass$gene_upper_without_dash, tfclass$gene_alter_name1_upper_without_dash,
                                                 tfclass$gene_alter_name2_upper_without_dash,tfclass$gene_alter_name3_upper_without_dash,
                                                 tfclass$gene_alter_name4_upper_without_dash,tfclass$gene_alter_name5_upper_without_dash,
                                                 tfclass$gene_alter_name6_upper_without_dash,tfclass$gene_alter_name7_upper_without_dash))


monomer_TFs_Lambert2018$Name_without_dash %in% unique_tfclass_gene_upper_without_dash %>% table()

#FALSE  TRUE
#71   647

# No correspondence found for these
# CENPB
# CPEB1
# SKOR2
# SKOR1
# ZNF704
# ZNF821
# MBNL2
# CPSF4
# DPF1
# XPA
# LIN28B
# ZNF385D 12


correspondence=read_delim("~/projects/TFBS/code/motif-networks/TFClass/Lambert2018_tfclass.txt", delim=" ")
correspondence[] <- lapply(correspondence, function(x) {
  if (is.character(x)) str_trim(x)
  else x
})

monomer_TF_Lambert2018_updated <- monomer_TFs_Lambert2018 %>%
  left_join(correspondence, by = c("Name_without_dash" = "Lambert2018")) %>%
  mutate(Name_without_dash = coalesce(TFClass, Name_without_dash)) %>%  # replace if mapping exists
  select(-TFClass)

monomer_TF_Lambert2018_updated$Name_without_dash %in% unique_tfclass_gene_upper_without_dash %>% table()

test=monomer_TF_Lambert2018_updated[which(!(monomer_TF_Lambert2018_updated$Name_without_dash %in% unique_tfclass_gene_upper_without_dash )),]
test$Name_without_dash #12

# Pivot table_b to long format (one gene name per row)
tfclass_long <- tfclass %>%
  pivot_longer(
    cols = ends_with("_without_dash"),
    names_to = "source_column",
    values_to = "gene_name"
  ) %>%
  filter(!is.na(gene_name))


monomer_TF_Lambert2018_updated=rename(monomer_TF_Lambert2018_updated, gene_name = Name_without_dash)

joined <- monomer_TF_Lambert2018_updated %>%
  left_join(tfclass_long, by = "gene_name")


which(table(joined$name)>1)
#MGA   MTF1   NRF1   TCF3   TCF4 ZSCAN1
#356    368    426    589    590    711

test=joined %>% filter(name %in% c("MGA",
                                   "MTF1",
                                   "NRF1",
                                   "TCF3",
                                   "TCF4",
                                   "ZSCAN1") )


df_combined <- joined %>%
  group_by(name) %>%
  summarise(
    superclass_code = toString(unique(superclass_code)),
    superclass = toString(unique(superclass)),
    class_code = toString(unique(class_code)),
    class = toString(unique(class)),
    class_abbr = toString(unique(class_abbr)),
    class_alter = toString(unique(class_alter)),
    family_code = toString(unique(family_code)),
    family.y = toString(unique(family.y)),
    family_abbr = toString(unique(family_abbr)),
    family_alter = toString(unique(family_alter)),
    subfamily_code = toString(unique(subfamily_code)),
    subfamily = toString(unique(subfamily)),
    subfamily_abbr = toString(unique(subfamily_abbr)),
    subfamily_alter = toString(unique(subfamily_alter)),
    gene_code = toString(unique(gene_code)),
    gene = toString(unique(gene)),
    gene_alter_name1 = toString(unique(gene_alter_name1)),
    gene_alter_name2 = toString(unique(gene_alter_name2)),
    gene_alter_name3 = toString(unique(gene_alter_name3)),
    gene_alter_name4 = toString(unique(gene_alter_name4)),
    gene_alter_name5 = toString(unique(gene_alter_name5)),
    gene_alter_name6 = toString(unique(gene_alter_name6)),
    gene_alter_name7 = toString(unique(gene_alter_name7)),
    numbers_after_gene = toString(unique(numbers_after_gene)),
    equals = toString(unique(equals)),
    numbers_after_equals = toString(unique(numbers_after_equals)),
    motif = toString(unique(motif)),
    uniprot = toString(unique(uniprot)),
    protein_atlas = toString(unique(protein_atlas)),
    pdb = toString(unique(pdb)),
    humanpsd = toString(unique(humanpsd)),
    transfac = toString(unique(transfac)),
    source_column = toString(unique(source_column)),

    .groups = "drop"
  )

df_combined$source_column=NULL


# Columns you want to merge
cols_to_merge <- c("gene_alter_name1", "gene_alter_name2", "gene_alter_name3",
                   "gene_alter_name4", "gene_alter_name5", "gene_alter_name6", "gene_alter_name7")


df_combined <- df_combined %>%
  mutate(gene_alter = pmap_chr(across(all_of(cols_to_merge)), ~ {
    vals <- c(...)               # collect values
    print(vals)
    vals <- vals[!is.na(vals) & vals != "NA"]  # remove NAs and empty strings
    if (length(vals) > 0) paste(vals, collapse = ", ") else NA_character_
  })) %>%
  relocate(gene_alter, .after = gene) %>%
  select(-all_of(cols_to_merge))  # remove original columns


cols_to_rename <- c(
  "superclass_code", "superclass", "class_code", "class", "class_abbr", "class_alter",
  "family_code", "family.y", "family_abbr", "family_alter",
  "subfamily_code", "subfamily", "subfamily_abbr", "subfamily_alter",
  "gene_code", "gene", "gene_alter", "numbers_after_gene",
  "equals", "numbers_after_equals", "motif", "uniprot", "protein_atlas",
  "pdb", "humanpsd", "transfac"
)

# Add "TFClass_" prefix to each column name
df_combined <- df_combined %>%
  rename_with(.fn = ~ paste0("TFClass_", .), .cols = all_of(cols_to_rename))

df_combined <- df_combined %>%
  rename(TFClass_family = TFClass_family.y)


#Combined with Lambert2018 data

monomer_TFs_Lambert2018_TFClass <- monomer_TF_Lambert2018_updated %>%
  left_join(df_combined, by = "name")

monomer_TFs_Lambert2018_TFClass$gene_name=NULL
monomer_TFs_Lambert2018_TFClass$Name_upper=NULL
monomer_TFs_Lambert2018_TFClass$Name_upper_without_dash=NULL

names(monomer_TFs_Lambert2018_TFClass)

lambert_cols <- c(
  "Name", "Ensembl_ID", "DBD", "Is TF?", "TF assessment", "Binding mode", "Motif status",
  "Notes", "Comments", "Committee notes", "MTW Notes", "TRH Notes", "SL notes", "AJ notes",
  "Disagree on Assessment", "Disagree on Binding",
  "Author1", "Assesment1", "Binding1", "Comment1", "Notes1",
  "Author2", "Assesment2", "Binding2", "Comment2", "Notes2",
  "Vaquerizas 2009 TF classification", "CisBP considers it as a TF?",
  "TFclass considers it as a TF?", "TF-CAT classification", "Is a GO TF"
)

# Apply prefix
monomer_TFs_Lambert2018_TFClass <- monomer_TFs_Lambert2018_TFClass%>%
  rename_with(~ paste0("Lambert2018_", .), .cols = all_of(lambert_cols))

#Add Wei 2010 ETS class

ETS_classes=read_csv("~/projects/cytoscape/monomer_network/ETS_classes.csv") %>% select(name, "shared name", "node_type", "ETS class")

ETS_classes= ETS_classes %>% filter(node_type=="TF" & !is.na(`ETS class`))

ETS_classes=ETS_classes %>% rename("Wei 2010 ETS class"="ETS class")


test<- monomer_TFs_Lambert2018_TFClass %>%
  left_join(ETS_classes %>% select(name, "Wei 2010 ETS class"), by = "name")

test=test %>% relocate(`Wei 2010 ETS class`, .after = TFClass_subfamily)

# Save the final data frame, these are only TF nodes

write_csv(test, "~/projects/cytoscape/monomer_network/node_table_monomer_TFs_Lambert2018_TFClass_ETS_class.csv")


