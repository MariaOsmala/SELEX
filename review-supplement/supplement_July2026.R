# region packages

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(httr2)
library(tibble)
#library(biomaRt)
#library(Biostrings)
library(httr)

# regionend


# region load data
rm(list = ls())

# All 3933 motifs with current type and protein copy info



metadata=read_delim("/Users/osmalama/projects/SELEX/motif-metadata/motifs_with_SELEX_signal_and_background_and_proteins_01072026.tsv",
            delim = "\t"
)



metadata |>
  filter(representative == "YES") |>
  dplyr::select(type) |>
  table(useNA = "always") |>
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

# metadata

data_path <- "/Users/osmalama/projects/SELEX/motif-metadata/"


metadata$TF1_copies |> table(useNA = "always")
metadata$TF2_copies |> table(useNA = "always")



# Check for which we need to create S2A data because they have protein info

is.na(metadata$`protein sequence 1`) |>
  table(useNA = "always") # 3326 with nonNA value

is.na(metadata$`TF1_protein_sequence`) |>
  table(useNA = "always") # 3608 with nonNA value


motifs_with_protein_info <- metadata |>
  filter(!is.na(`protein sequence 1`)) |>
  dplyr::select(ID) #3326

write_delim(motifs_with_protein_info,
  paste0(
    data_path,
    "motifs_with_protein_info_01072026.tsv"
  ),
  delim = "\t"
)

# motifs with S2A data, those with protein_info and k-mer scores? no, less

# metadata %>% filter(!is.na(`protein sequence 1`) & kmer_scoring_exists==TRUE) %>% nrow() #3290 matches

# included_in_S2A=metadata %>% filter(!is.na(`protein sequence 1`) & kmer_scoring_exists==TRUE) %>% pull(ID)

# metadata$included_in_S2A=FALSE

# metadata$included_in_S2A[which(metadata$ID %in% included_in_S2A)]=TRUE

#use this once ready
#"/projappl/project_2013895/motif_metadata/motifs_with_S2A_and_proteins_01072026.tsv"

data_file <- "motifs_with_S2A_and_proteins_01072026.tsv"
s2a <- read_delim(paste0(data_path, data_file), delim = "\t") #

s2a$included_in_S2A |> table(useNA = "always") # 
# FALSE  TRUE  <NA>
#  656  3277     0


metadata <- metadata |>
  left_join(
    s2a |>
      dplyr::select(ID, included_in_S2A),
    by = "ID"
  )

metadata$included_in_S2A |> table(useNA = "always") # 3277
#FALSE  TRUE  <NA> 
#  656  3277     0 

#What are those for which k-mer scores and proteins exist but do not end up to the final dataset

metadata %>% filter(!is.na(`protein sequence 1`) & kmer_scoring_exists==TRUE & !included_in_S2A) %>% select(ID, TF1_copies, TF2_copies)

# Some issue with this data
# 1 MGA_DLX2_CAP-SELEX_TGCGGT40NTCA_AAD_AGGTGNTAATTR_1_2u             1          1         
# 2 POU2F1_ELK1_CAP-SELEX_TGACGA40NGCA_AS_NCCGGATATGCAN_1_2u          1          1         
# 3 POU2F1_ETV1_CAP-SELEX_TGCGAA40NAGC_AS_NCCGGATATGCAN_1_2u          1          1         
# 4 ERF_SREBF2_CAP-SELEX_TCTTTG40NACT_AAC_NNCACGTGACMGGAARNN_1_3u     1          2         
# 5 GCM1_ERG_CAP-SELEX_TAAGAA40NATA_AX_RTRYGGGCGGAARKN_1_3u           1          1         
# 6 HOXA3_PAX5_CAP-SELEX_TCTGTC40NCAA_AY_YNATTAGTCACGCWTSRNTR_2_3u    1          1         
# 7 GCM2_DLX2_CAP-SELEX_TGTGCT40NCGG_AAB_RTRCGGGNNNNNTAATTR_1_3u      1          1         
# 8 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_ATRCGGGYNNNNNYWTTGTNN_1_3u  1          1         
# 9 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNNNNNYWTTGTNN_1_3u 1          1         
# 10 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNNRNACAAWN_1_3u    1          1         
# 11 GCM2_SOX15_CAP-SELEX_TCGCCA40NCCT_AAB_RTRCGGGNNNRNACAAWN_1_3u     1          1         
# 12 HOXB7_HT-SELEX_TTATTT40NCGT_KV_NGTAATTANN_1_4u                    1          NA        
# 13 TBX18_HT-SELEX_TATTTG40NCCA_KP_NRAGGTGTGAAN_1_4u                  1          NA     

# S2A download link
# https://a3s.fi/SELEX/sequence-to-affinity-data.tar.gz


metadata |>
  dplyr::select(type) |>
  table(useNA = "always") # 

metadata |>
  filter(representative == "YES") |>
  dplyr::select(type) |>
  table(useNA = "always") # 

metadata %>% filter(experiment!="CAP-SELEX" & representative=="YES") %>% select(TF1_copies) %>% table(useNA="always") 
#  1    2    3    4 <NA> 
#  230  186    3    9   11 

# Convert YES/NO in representative to yes/no

metadata <- metadata |>
  mutate(representative = ifelse(representative == "YES", "yes",
    ifelse(representative == "NO", "no", representative)
  ))

# Convert TRUE FALSE to yes/no

metadata <- metadata |>
  mutate(included_in_S2A = ifelse(included_in_S2A == TRUE, "yes",
    ifelse(included_in_S2A == FALSE, "no", included_in_S2A)
  ))


metadata |>
  filter(representative == "yes") |>
  dplyr::select(type) |>
  table(useNA = "always") %>% sum()# 

metadata<-  metadata%>%
  separate(Lambert2018_families, into = c("TF1_family", "TF2_family"), sep = "_")

# TF1_family
# TF2_family

# Relocate TF1_family and TF2_family after TF1 and TF2

metadata <- metadata %>%
  relocate(TF1_family, .after = TF1) %>%
  relocate(TF2_family, .after = TF2)

metadata %>% filter(experiment!="CAP-SELEX"& is.na(TF1_family)) %>% pull(ID)
#THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3 Homeodomain

metadata$TF1_family[which(metadata$ID=="THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3")]="Homeodomain"

metadata %>% filter(experiment!="CAP-SELEX"& TF1_family=="Znf") %>% pull(ID)
#"XPA_HT-SELEX_TCACGC40NACT_KX_NCACCTCACAN_1_4" Znf_XPA
metadata$TF1_family[which(metadata$ID=="XPA_HT-SELEX_TCACGC40NACT_KX_NCACCTCACAN_1_4")]="Znf_XPA"


#write_delim(metadata, file="/Users/osmalama/projects/TFBS/Data/SELEX-motif-collection/Supplementary_Table_1_Submission_March2026.tsv", delim="\t")
write_delim(metadata, file="/Users/osmalama/projects/TFBS/Data/SELEX-motif-collection/Supplementary_Table_1_Submission_July2026.tsv", delim="\t")
#Move this to "/Users/osmalama/Dropbox/Taipale-lab/Nature genetics review/Submission_July2026/metadata.tsv"
# Write table to excel

# install.packages("openxlsx")
library(openxlsx)


# --- Export in the "Base + A/C/G/T" block format ---
# We need TF1_protein_sequence, TF2_protein_sequence


pfm_path_col <- "filename"
out_file <- "/Users/osmalama/Dropbox/Taipale-lab/Nature genetics review/Submission_July2026/Supplementary_Table_1.xlsx" # nolint
export_cols <- c(
  "motif_ID", "experiment", "study", "TF1", "TF2",
  "TF1_family", "TF2_family", "type", "representative",
  "TF1_protein_sequence", "TF2_protein_sequence",
  "TF1_copies", "TF2_copies", "seed",
  "included_in_S2A","fastq_ftp_signal", "fastq_ftp_background"   
)

            

add_blank_row_between <- FALSE
fontName <- "Verdana"
fontSize <- 12

# Read all PFMs
setwd("/Users/osmalama/projects/TFBS/RProjects/TFBS") # nolint
pfms <- lapply(metadata[[pfm_path_col]],
  read_delim,
  delim = "\t", col_names = FALSE
)
motif_lengths <- vapply(pfms, ncol, integer(1))
maxL <- max(motif_lengths, na.rm = TRUE)

# If motifs are longer than metadata column count, add extra position columns
extra_needed <- max(0, maxL - length(export_cols))
pos_cols <- if (extra_needed > 0) {
  sprintf("pos_%02d", (length(export_cols) + 1):maxL)
} else {
  character(0)
}

out_cols <- c(export_cols, pos_cols)


# Workbook
wb <- createWorkbook()
addWorksheet(wb, "motifs")

# Header row (blue)
header <- c("", out_cols)

writeData(wb, "motifs",
  x = as.data.frame(t(header)),
  startRow = 1, startCol = 1, colNames = FALSE
)

headerStyle <- createStyle( # nolint
  fontName = fontName, fontSize = fontSize,
  fgFill = "#00B0F0", textDecoration = "bold",
  halign = "center", valign = "center"
)
baseStyle <- createStyle( # nolint
  fontName = fontName, fontSize = fontSize,
  fgFill = "#92D050"
)

fontStyle <- createStyle(fontName = fontName, fontSize = fontSize) # nolint

addStyle(wb, "motifs", headerStyle,
  rows = 1, cols = 1:length(header),
  gridExpand = TRUE
)

r <- 2
for (i in seq_len(nrow(metadata))) {
  # --- Base row label (col 1) ---
  writeData(wb, "motifs", "Base", startRow = r, startCol = 1, colNames = FALSE)

  # --- Base row metadata (cols 2..): keep original types ---
  base_row <- metadata[i, export_cols, drop = FALSE]
  # add extra position columns as NA_real_ so they are numeric in Excel (blank)
  for (pc in pos_cols) base_row[[pc]] <- NA_real_
  base_row <- base_row[, out_cols, drop = FALSE]

  writeData(wb, "motifs", base_row,
    startRow = r, startCol = 2, colNames = FALSE
  )

  # style Base row green
  addStyle(wb, "motifs", baseStyle,
    rows = r,
    cols = 1:(length(out_cols) + 1), gridExpand = TRUE, stack = TRUE
  )

  # --- A/C/G/T labels (col 1) ---
  writeData(wb, "motifs", c("A", "C", "G", "T"),
    startRow = r + 1, startCol = 1, colNames = FALSE
  )

  # --- Counts numeric matrix (cols 2..): ALWAYS numeric in Excel ---
  pfm <- as.matrix(pfms[[i]])
  rownames(pfm) <- c("A", "C", "G", "T")
  L <- ncol(pfm)

  counts <- matrix(NA_real_, nrow = 4, ncol = length(out_cols))
  rownames(counts) <- c("A", "C", "G", "T")
  counts[, 1:L] <- pfm[rownames(counts), , drop = FALSE]

  writeData(wb, "motifs", counts,
    startRow = r + 1, startCol = 2,
    colNames = FALSE, rowNames = FALSE
  )

  # advance row pointer
  r <- r + 5
  if (add_blank_row_between) r <- r + 1
}

# Apply Verdana 12 everywhere we wrote (stack so it doesn't wipe fills)
addStyle(wb, "motifs", fontStyle,
  rows = 1:(r - 1), cols = 1:(length(out_cols) + 1),
  gridExpand = TRUE, stack = TRUE
)

freezePane(wb, "motifs", firstRow = TRUE)
# setColWidths(wb, "motifs", cols = 1:(length(out_cols) + 1), widths = "auto")

saveWorkbook(wb, out_file, overwrite = TRUE)
out_file

#What is the stuff below? Results to exactly the same excel table?

# Build output table
blocks <- vector("list", nrow(metadata))

for (i in seq_len(nrow(metadata))) {
  meta <- metadata[i, export_cols, drop = FALSE]
  pfm <- as.matrix(pfms[[i]])
  rownames(pfm) <- c("A", "C", "G", "T")
  L <- ncol(pfm)

  # Base row = metadata
  base_row <- as.list(c(row_type = "Base", meta))
  for (pc in pos_cols) base_row[[pc]] <- NA

  # A/C/G/T rows = counts filled left-to-right across out_cols
  acgt_rows <- lapply(c("A", "C", "G", "T"), function(b) {
    vals <- rep(NA_real_, length(out_cols))
    vals[seq_len(L)] <- as.numeric(pfm[b, ])
    c(list(row_type = b), as.list(setNames(vals, out_cols)))
  })

  block_metadata <- do.call(rbind, c(
    list(as.data.frame(base_row, stringsAsFactors = FALSE)),
    lapply(acgt_rows, as.data.frame, stringsAsFactors = FALSE)
  ))

  if (add_blank_row_between) {
    blank <- as.data.frame(
      c(
        list(row_type = NA_character_),
        as.list(setNames(rep(NA, length(out_cols)), out_cols))
      ),
      stringsAsFactors = FALSE
    )
    block_metadata <- rbind(block_metadata, blank)
  }
  blocks[[i]] <- block_metadata
}

out_df <- do.call(rbind, blocks)

# Write Excel + styles
wb <- createWorkbook()
addWorksheet(wb, "motifs")

# Make first header cell blank (like your screenshot)
colnames(out_df)[1] <- ""

writeData(wb, "motifs", out_df, colNames = TRUE)

n_rows <- nrow(out_df) + 1 # +1 for header row
n_cols <- ncol(out_df)

# Default font for everything
fontStyle <- createStyle(fontName = "Verdana", fontSize = 12)
addStyle(wb, "motifs", fontStyle,
  rows = 1:n_rows, cols = 1:n_cols,
  gridExpand = TRUE, stack = TRUE
)

# If you also want header/base rows to keep their fills AND use Verdana 12:
headerStyle <- createStyle(
  fontName = "Verdana", fontSize = 12,
  fgFill = "#00B0F0", textDecoration = "bold",
  halign = "center", valign = "center"
)
baseStyle <- createStyle(
  fontName = "Verdana", fontSize = 12,
  fgFill = "#92D050"
)

addStyle(wb, "motifs", headerStyle,
  rows = 1, cols = 1:n_cols,
  gridExpand = TRUE, stack = TRUE
)

base_rows <- which(out_df[[1]] == "Base") + 1
addStyle(wb, "motifs", baseStyle,
  rows = base_rows, cols = 1:n_cols,
  gridExpand = TRUE, stack = TRUE
)


base_rows <- which(out_df[[1]] == "Base") + 1 # +1 for header row
if (length(base_rows)) {
  addStyle(wb, "motifs", baseStyle,
    rows = base_rows, cols = 1:ncol(out_df),
    gridExpand = TRUE, stack = TRUE
  )
}

freezePane(wb, "motifs", firstRow = TRUE)
setColWidths(wb, "motifs", cols = 1:ncol(out_df), widths = "auto")

saveWorkbook(wb,"/Users/osmalama/Dropbox/Taipale-lab/Nature genetics review/Submission_July2026/Supplementary_Table_1_version2.xlsx", overwrite = TRUE)
out_file


# ------------------ Example usage ------------------
# Suppose your data.frame is called `metadata` and has a column with .pfm paths called "pfm_path"
# out <- export_motifs_to_xlsx(metadata, pfm_path_col = "pfm_path", out_file = "my_motifs.xlsx")
# out
