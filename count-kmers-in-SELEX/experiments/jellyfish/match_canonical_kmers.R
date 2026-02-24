#!/usr/bin/env Rscript

reverse_complement <- function(seq) {
  s <- paste(rev(strsplit(seq, "", fixed = TRUE)[[1]]), collapse = "")
  chartr("ACGTN", "TGCAN", s)
}

read_lines_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    readLines(con)
  } else {
    readLines(path)
  }
}

load_canonical_counts <- function(counts_file) {
  lines <- read_lines_any(counts_file)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\\s+")
  ok <- vapply(parts, length, integer(1)) == 2L
  if (!all(ok)) parts <- parts[ok]
  
  kmers  <- vapply(parts, `[`, character(1), 1)
  counts <- as.integer(vapply(parts, `[`, character(1), 2))
  
  # Aggregate duplicates (shouldn’t exist in jellyfish dump, but safe)
  df <- data.frame(kmer = kmers, count = counts, stringsAsFactors = FALSE)
  agg <- rowsum(df$count, df$kmer)
  out <- as.integer(agg[, 1])
  names(out) <- rownames(agg)
  out
}

match_pattern_canonical <- function(pattern, counts_vec) {
  rc_pattern <- reverse_complement(pattern)
  forward_regex <- paste0("^", gsub("N", ".", pattern, fixed = TRUE), "$")
  reverse_regex <- paste0("^", gsub("N", ".", rc_pattern, fixed = TRUE), "$")
  
  kmers <- names(counts_vec)
  m_fwd <- grepl(forward_regex, kmers, perl = TRUE)
  m_rev <- grepl(reverse_regex, kmers, perl = TRUE)
  idx   <- which(m_fwd | m_rev)
  
  total_count  <- if (length(idx)) sum(counts_vec[idx]) else 0L
  matched_kmers <- if (length(idx)) kmers[idx] else character(0)
  
  if (length(matched_kmers)) {
    rc_matched <- vapply(matched_kmers, reverse_complement, character(1))
    matched_kmers <- unique(c(matched_kmers, rc_matched))
  }
  
  list(total_count = as.integer(total_count), matched_kmers = matched_kmers)
}

process_files <- function(kmers_file, background_counts, signal_counts, output_file) {
  
  # kmers_file="/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3.txt"
  # background_counts="/scratch/project_2013895/SELEX/jellyfish/filtered_counts_all_kmers/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3_background_counts.txt"
  # signal_counts="/scratch/project_2013895/SELEX/jellyfish/filtered_counts_all_kmers/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3_signal_counts.txt"
  # output_file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/jellyfish/output.txt"
  
  # Load k-mers with Ns
  patterns <- trimws(read_lines_any(kmers_file))
  patterns <- patterns[nzchar(patterns)]
  
  cat("Loading background counts...\n")
  bg_counts <- load_canonical_counts(background_counts)
  
  cat("Loading signal counts...\n")
  sig_counts <- load_canonical_counts(signal_counts)
  
  cat(sprintf("Processing %d patterns...\n", length(patterns)))
  
  results_list <- vector("list", length(patterns))
  fold_vals <- numeric(length(patterns))
  
  for (i in seq_along(patterns)) {
    p <- patterns[i]
    
    bg <- match_pattern_canonical(p, bg_counts)
    sg <- match_pattern_canonical(p, sig_counts)
    
    bg_count <- bg$total_count
    sg_count <- sg$total_count
    
    fold_enrichment <- if (bg_count > 0L) {
      sg_count / bg_count
    } else if (sg_count > 0L) {
      Inf
    } else {
      NA_real_
    }
    fold_vals[i] <- fold_enrichment
    
    results_list[[i]] <- list(
      pattern          = p,
      rc_pattern       = reverse_complement(p),
      background_count = bg_count,
      signal_count     = sg_count,
      fold_enrichment  = fold_enrichment,
      n_matched_kmers  = length(union(bg$matched_kmers, sg$matched_kmers))
    )
  }
  
  # Prepare output table (format fold to 4 decimals when finite)
  fmt_fold <- function(x) {
    if (is.na(x)) return("NA")
    if (is.infinite(x)) return("Inf")
    formatC(x, format = "f", digits = 4)
  }
  
  out <- data.frame(
    `K-mer`             = vapply(results_list, `[[`, character(1), "pattern"),
    RC_Pattern          = vapply(results_list, `[[`, character(1), "rc_pattern"),
    Background_Count    = vapply(results_list, `[[`, integer(1),   "background_count"),
    Signal_Count        = vapply(results_list, `[[`, integer(1),   "signal_count"),
    Fold_Enrichment     = vapply(fold_vals, fmt_fold, character(1)),
    N_Matched_Kmers     = vapply(results_list, `[[`, integer(1),   "n_matched_kmers"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  write.table(out, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Results saved to %s\n", output_file))
  
  # Summary
  enriched <- sum(is.finite(fold_vals) & !is.na(fold_vals) & fold_vals > 1)
  depleted <- sum(is.finite(fold_vals) & !is.na(fold_vals) & fold_vals < 1)
  
  cat("\nSummary:\n")
  cat(sprintf("Total patterns processed: %d\n", length(patterns)))
  cat(sprintf("Enriched in signal: %d\n", enriched))
  cat(sprintf("Depleted in signal: %d\n", depleted))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4L) {
    #Rscript match_canonical_kmers.R kmers_with_N.txt background_counts.txt signal_counts.txt output.txt
    cat("Usage: Rscript match_canonical_kmers.R <kmers_with_N.txt> <background_counts.txt> <signal_counts.txt> <output.txt>\n")
    quit(status = 1)
  }
  # args=c( "/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3.txt",
  #          "/scratch/project_2013895/SELEX/jellyfish/filtered_counts_all_kmers/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3_background_counts.txt",
  #          "/scratch/project_2013895/SELEX/jellyfish/filtered_counts_all_kmers/NFATC4_ONECUT1_TCCGTC40NGTT_YJII_NRYGGAAANNNATCGATYKN_2_3_signal_counts.txt",
  #          "/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/jellyfish/output.txt")
  process_files(args[1], args[2], args[3], args[4])
}

if (identical(environment(), globalenv())) {
  main()
}
