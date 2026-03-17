# region packages

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(httr2)
library(tibble)
library(biomaRt)
library(Biostrings)
library(httr)

# regionend


# region load data
rm(list = ls())

# All 3933 motifs with current type and protein copy info

# Clean metadata

metadata <- read_delim("/Users/osmalama/projects/TFBS/Data/SELEX-motif-collection/metadata_final.tsv", delim = "\t") # nolint #3933

included_cols <- c(
  "ID", "motif_ID", "symbol", "study", "experiment", "seed", "multinomial", # nolint
  "cycle", "representative", "filename", "IC", "length", # nolint
  "Methyl.SELEX.Motif.Category" # nolint
) # nolint


metadata <- metadata |> dplyr::select(all_of(included_cols))


# Monomer network

monomer_nodes <- read_delim(
  "~/projects/cytoscape/monomer_network/node_table_24022026.csv", # nolint
  col_types = cols(.default = col_character()), delim = ","
) # nolint

monomer_nodes <- monomer_nodes |> filter(node_type == "motif") # 2035

monomer_nodes |>
  filter(representative == "YES") |>
  nrow() # 439

included_cols <- c(
  "name", "symbol",
  "Lambert2018_families",
  "type",
  "type_review"
)

monomer_nodes <- monomer_nodes |>
  dplyr::select(all_of(included_cols)) |>
  dplyr::rename(ID = name, TF1 = symbol, TF1_family = Lambert2018_families)

monomer_nodes$TF2 <- NA
monomer_nodes$TF2_family <- NA


# Heterodimer network

heterodimer_nodes <- read_delim(
  "~/projects/cytoscape/heterodimer_network/node_table_24022026.csv", # nolint
  col_types = cols(.default = col_character()), delim = ","
) # nolint

heterodimer_nodes <- heterodimer_nodes |> filter(node_type == "motif") # 1898

heterodimer_nodes |>
  filter(representative == "YES") |>
  nrow() # 793


included_cols <- c(
  "name",
  "TF1", "TF2", "TF1_families",
  "TF2_families", "type", "type_review"
)


heterodimer_nodes <- heterodimer_nodes |>
  dplyr::select(all_of(included_cols)) |>
  dplyr::rename(ID = name, TF1_family = TF1_families, TF2_family = TF2_families)

node_table <- rbind(
  monomer_nodes |>
    dplyr::select(ID, TF1, TF1_family, TF2, TF2_family, type, type_review),
  heterodimer_nodes |>
    dplyr::select(ID, TF1, TF1_family, TF2, TF2_family, type, type_review)
)

metadata <- metadata |> left_join(node_table, by = "ID")

metadata$type[which(is.na(metadata$type))] <- metadata$type_review[which(is.na(metadata$type))] # nolint

metadata |>
  filter(representative == "YES") |>
  dplyr::select(type) |>
  table(useNA = "always") |>
  as.data.frame()

#              type Freq
# 1       composite  483
# 2      composite?    1
# 3         dimeric  179
# 4        dimeric?    9
# 5       monomeric  211
# 6         spacing  301
# 7        spacing?    8
# 8   ssDNA binding   17
# 9  ssDNA binding?    2
# 10     tetrameric    6
# 11    tetrameric?    3
# 12       trimeric    3
# 13        unknown    9
# 14           <NA>    0

# metadata

metadata |>
  filter(representative == "YES") |>
  dplyr::select(type_review) |>
  table(useNA = "always") |>
  as.data.frame()

# type_review Freq
# 1     composite  489
# 2       dimeric  189
# 3     monomeric  211
# 4       spacing  304
# 5 ssDNA binding   18
# 6    tetrameric    8
# 7      trimeric    4
# 8       unknown    9
# 9          <NA>    0

data_path <- "/Users/osmalama/projects/SELEX/motif-metadata/"

# motifs for which there is SELEX data and the protein sequence is checked against the claimed ENSG
data_file <-
  "motifs_with_SELEX_signal_and_background_and_proteins_04032026_checked.tsv"

proteins_checked <- read_delim(paste0(data_path, data_file),
  col_types = cols(.default = col_character()),
  delim = "\t"
) # 3224

proteins_checked$TF1_copies |> table(useNA = "always")
proteins_checked$TF2_copies |> table(useNA = "always")


# test that TF1 and TF2 are the same

match_ind <- match(proteins_checked$ID, metadata$ID)

table(proteins_checked$TF1 == metadata$TF1[match_ind]) # TRUE
table(proteins_checked$TF2 == metadata$TF2[match_ind]) # TRUE

metadata <- metadata |>
  left_join(
    proteins_checked |>
      dplyr::select( # nolint
        ID, TF1_copies, TF2_copies, # nolint
        TF1_protein_sequence, TF2_protein_sequence,
        `protein sequence 1`, `protein sequence 2`, # nolint
        `protein sequence 3`, `protein sequence 4`
      ), # nolint
    by = "ID" # nolint
  )

# Check for which we need to create S2A data because they have protein info

is.na(metadata$`protein sequence 1`) |>
  table(useNA = "always") # 3224 with nonNA value

is.na(metadata$`TF1_protein_sequence`) |>
  table(useNA = "always") # 3224 with nonNA value


motifs_with_protein_info <- metadata |>
  filter(!is.na(`protein sequence 1`)) |>
  dplyr::select(ID)

write_delim(motifs_with_protein_info,
  paste0(
    data_path,
    "motifs_with_protein_info_04032026.tsv"
  ),
  delim = "\t"
)

# motifs with S2A data

data_file <- "motifs_with_S2A_and_proteins.tsv"
s2a <- read_delim(paste0(data_path, data_file), delim = "\t") # 3594

s2a$included_in_S2A |> table(useNA = "always") # 3594 with nonNA value
# FALSE  TRUE  <NA>
#  487  3107     0


metadata <- metadata |>
  left_join(
    s2a |>
      dplyr::select(ID, included_in_S2A),
    by = "ID"
  )

metadata$included_in_S2A |> table(useNA = "always") # 3107
# FALSE  TRUE  <NA>
#  487  3107   339

metadata$included_in_S2A[which(is.na(metadata$included_in_S2A))] <- FALSE
# FALSE  TRUE  <NA>
# 826  3107     0

# S2A download link
# https://a3s.fi/SELEX/sequence-to-affinity-data.tar.gz


metadata |>
  dplyr::select(type_review) |>
  table(useNA = "always") # 483 composite, 179 dimeric, 211 monomeric, 301 spacing, 17 ssDNA binding, 6 tetrameric, 3 trimeric, 9 unknown, 8 NA

metadata |>
  filter(representative == "YES") |>
  dplyr::select(type_review) |>
  table(useNA = "always") # 483 composite, 179 dimeric, 211 monomeric, 301 spacing, 17 ssDNA binding, 6 tetrameric, 3 trimeric, 9 unknown, 8 NA


# For 2015 CAP-SELEX non-representative motifs, convert the type_review to NA

metadata <- metadata |>
  mutate(type_review = ifelse((experiment == "CAP-SELEX" & study == "Jolma2015" & representative == "NO"), NA, type_review))

metadata$type <- NULL

# rename type_review to type

metadata <- metadata |> dplyr::rename("type" = "type_review")

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


# Write table to excel

# install.packages("openxlsx")
library(openxlsx)


# --- Export in the "Base + A/C/G/T" block format ---
# We need TF1_protein_sequence, TF2_protein_sequence


pfm_path_col <- "filename"
out_file <- "/Users/osmalama/Dropbox/Taipale-lab/Nature genetics review/Submission_March2026/Supplementary_Table_1.xlsx" # nolint
export_cols <- c(
  "motif_ID", "experiment", "study", "TF1", "TF2",
  "TF1_family", "TF2_family", "type", "representative",
  "TF1_protein_sequence", "TF2_protein_sequence",
  "TF1_copies", "TF2_copies", "seed",
  "included_in_S2A"
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

saveWorkbook(wb, out_file, overwrite = TRUE)
out_file


# ------------------ Example usage ------------------
# Suppose your data.frame is called `metadata` and has a column with .pfm paths called "pfm_path"
# out <- export_motifs_to_xlsx(metadata, pfm_path_col = "pfm_path", out_file = "my_motifs.xlsx")
# out
