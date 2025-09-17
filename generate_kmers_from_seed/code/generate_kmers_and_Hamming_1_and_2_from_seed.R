#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(Biostrings)
  library(TFBSTools)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
arrays <- as.numeric(args[1])
addition=as.numeric(args[2])
arrays=arrays+addition



## ------------------------------------------------------------------
## Inputs/paths
## ------------------------------------------------------------------
metadata <- read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")
metadata <- metadata %>% filter(!(experiment %in% "Methyl-HT-SELEX"))
results_path <- "/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/"
pseudocount <- 0.01

## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
canonical_form <- function(seq) {
  s  <- DNAString(seq)
  rc <- as.character(reverseComplement(s))
  s  <- as.character(s)
  if (s <= rc) s else rc
}

## Generate Hamming neighbourhood up to radius 2 (substitutions only)
generate_hamming <- function(kmer, radius = 2) {
  chars <- strsplit(kmer, "", fixed = TRUE)[[1]]
  n <- length(chars)
  bases <- c("A","C","G","T")
  out <- character(0)
  
  ## radius 1
  if (radius >= 1) {
    for (i in seq_len(n)) {
      for (alt in bases[bases != chars[i]]) {
        m1 <- chars; m1[i] <- alt
        out <- c(out, canonical_form(paste0(m1, collapse = "")))
      }
    }
  }
  
  ## radius 2
  if (radius >= 2 && n >= 2) {
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        for (alt1 in bases[bases != chars[i]]) {
          for (alt2 in bases[bases != chars[j]]) {
            m2 <- chars
            m2[i] <- alt1; m2[j] <- alt2
            out <- c(out, canonical_form(paste0(m2, collapse = "")))
          }
        }
      }
    }
  }
  
  unique(out)
}

## faster matrix-index scoring
score_sequence_fast <- function(seq, pssm_mat, base_index) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  rows  <- base_index[chars]
  sum(pssm_mat[cbind(rows, seq_along(rows))])
}

## ------------------------------------------------------------------
## Stream canonical k-mers + Hamming-1 + Hamming-2 neighbors (deduped)
## ------------------------------------------------------------------
stream_canonical_hamm12_to_file <- function(iupac_seq, output_file) {
  iupac_map <- Biostrings::IUPAC_CODE_MAP
  chars <- unlist(strsplit(iupac_seq, "", fixed = TRUE))
  base_options <- lapply(chars, function(ch) strsplit(iupac_map[[ch]], "")[[1]])
  
  seen <- new.env(hash = TRUE, parent = emptyenv())
  con <- file(output_file, open = "w")
  buffer <- character(0L)
  buf_cap <- 10000L
  flush_buffer <- function() {
    if (length(buffer)) {
      writeLines(buffer, con)
      buffer <<- character(0L)
    }
  }
  on.exit({ flush_buffer(); close(con) }, add = TRUE)
  
  stream_recursive <- function(pos = 1L, prefix = "") {
    if (pos > length(base_options)) {
      kmer  <- prefix
      canon <- canonical_form(kmer)
      
      if (!exists(canon, envir = seen, inherits = FALSE)) {
        assign(canon, TRUE, envir = seen)
        buffer <<- c(buffer, canon)
        if (length(buffer) >= buf_cap) flush_buffer()
      }
      
      vset <- generate_hamming(kmer, radius = 2)  # includes radius 1 & 2
      for (v in vset) {
        if (!exists(v, envir = seen, inherits = FALSE)) {
          assign(v, TRUE, envir = seen)
          buffer <<- c(buffer, v)
          if (length(buffer) >= buf_cap) flush_buffer()
        }
      }
    } else {
      for (b in base_options[[pos]]) {
        stream_recursive(pos + 1L, paste0(prefix, b))
      }
    }
  }
  
  stream_recursive()
}

## ------------------------------------------------------------------
## Main
## ------------------------------------------------------------------

data <- metadata[arrays, ]


## Scratch paths
path_local_scratch <- Sys.getenv("LOCAL_SCRATCH")


kmer_file  <- file.path(path_local_scratch, paste0(data$ID, ".tsv"))
score_file <- file.path(path_local_scratch, paste0(data$ID, "_score.tsv"))
rank_file  <- file.path(path_local_scratch, paste0(data$ID, "_score_rank.tsv"))

## 1) Stream k-mers (canonical + Hamming-1 + Hamming-2)
message("Streaming k-mers (seed + H1 + H2) for seed: ", data$seed)
stream_canonical_hamm12_to_file(data$seed, kmer_file)

#C is PFM, p is pseudocount
PPMp <- function(C, p){
  (C + p / nrow(C)) / matrix( (colSums(C) + p),nrow=4, ncol=ncol(C),byrow=TRUE)
  
}

# Position-specific scoring matrix ----------------------------------------

#Position specific scoring matrix
PSSM <- function(C, B) log(C / B) #ln


## 2) Build PSSM
pcm <- as.matrix(read.table(file = gsub("../../","/projappl/project_2006203/TFBS/", data$filename)))
dimnames(pcm) <- list(c("A","C","G","T"))
pfm <- TFBSTools::PFMatrix(
  ID = data$ID,
  strand = "+",
  bg = c(A=.25, C=.25, G=.25, T=.25),
  profileMatrix = pcm
)


PWM=PPMp(pfm@profileMatrix,pseudocount)

pssm=PSSM(PWM, 0.25)

base_index <- c(A=1L, C=2L, G=3L, T=4L)

## 3) Score k-mers by motif (chunked + optional parallel)
score_kmers_with_motif <- function(kmer_file, pssm, output_file, chunk = 1e5, mc.cores = NULL) {
  if (is.null(mc.cores)) mc.cores <- max(1L, detectCores() - 1L)
  in_con  <- file(kmer_file, open = "r")
  out_con <- file(output_file, open = "w")
  on.exit({ close(in_con); close(out_con) }, add = TRUE)
  
  repeat {
    lines <- readLines(in_con, n = chunk, warn = FALSE)
    if (!length(lines)) break
    
    if (mc.cores > 1L) {
      scored <- mclapply(lines, function(line) {
        rc <- as.character(reverseComplement(DNAString(line)))
        best <- max(
          score_sequence_fast(line, pssm, base_index),
          score_sequence_fast(rc,   pssm, base_index)
        )
        paste(line, best, sep = "\t")
      }, mc.cores = mc.cores)
      writeLines(unlist(scored, use.names = FALSE), out_con)
    } else {
      for (line in lines) {
        rc <- as.character(reverseComplement(DNAString(line)))
        best <- max(
          score_sequence_fast(line, pssm, base_index),
          score_sequence_fast(rc,   pssm, base_index)
        )
        writeLines(paste(line, best, sep = "\t"), out_con)
      }
    }
  }
}

message("Scoring k-mers …")
score_kmers_with_motif(kmer_file = kmer_file, pssm = pssm, output_file = score_file)

system(paste("rm", shQuote(kmer_file)))

#result=read_tsv(paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv"))
#result=read_tsv(paste0(path_local_scratch,"/",metadata$ID[arrays],"_score_rank.tsv"), col_names=NULL)
result_old=read_tsv(
paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX","/",metadata$ID[arrays],".tsv"),
col_names=NULL)

## 4) Sort by score desc, add rank, and write final artifact
cmd <- sprintf(
  "sort -k2,2nr %s | awk '{printf \"%%s\\t%%s\\t%%d\\n\", $1, $2, NR}' > %s",
  shQuote(score_file), shQuote(rank_file)
)
system(cmd)
system(paste("rm", shQuote(score_file)))

## 5) Copy to results path
system(paste("cp", shQuote(rank_file), shQuote(results_path)))

message("Done. Final ranked file: ", rank_file)
