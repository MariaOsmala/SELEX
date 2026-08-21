#!/usr/bin/env Rscript

# Suppress the Biostrings startup chatter so the script output stays clean in pipelines:
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse))

# Rscript code/lambda.R --sig ../data/FOXO1_TTCAGC20NTA_AF_3.seq \
# --bg  ../data/FOXO1_TTCAGC20NTA_AF_2.seq \
# --length 20 --kmer 8 --out lambda_output.tsv > test.tsv
# 
# Rscript lambda.R --help     # auto-generated usage message


count_kmer_middle <- function(path, k = 8, expected_length = NULL) {
  # --- 1. Read sequences based on file extension ---
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("fa", "fasta", "fna", "ffn")) {
    seqs_raw <- readDNAStringSet(path, format = "fasta")
  } else if (ext %in% c("fq", "fastq", "gz")) {
    seqs_raw <- readDNAStringSet(path, format = "fastq")
  } else {
    # .seq or unknown: one sequence per line
    lines <- toupper(trimws(readLines(path))) # Strips leading and trailing whitespace from each string. Converts everything to uppercase
    lines <- lines[nzchar(lines)] # returns TRUE for non-empty strings (
    seqs_raw <- DNAStringSet(lines)
  }
  
  n_orig <- length(seqs_raw) #225057
  
  # --- 2. Length filter (applied to all formats) ---
  if (!is.null(expected_length)) {
    seqs_raw <- seqs_raw[width(seqs_raw) == expected_length]
  }
  
  n_after_length_filtering <- length(seqs_raw)
  
  # --- 2. N / ambiguity filtering (reject anything not pure ACGT) ---
  amb_counts <- letterFrequency(seqs_raw, letters = "NRYSWKMBDHV")
  seqs <- seqs_raw[rowSums(amb_counts) == 0]
  n_after_ambiguity_filtering <- length(seqs) #225026
  
  # --- 3. Add reverse complements ---
  seqs2 <- c(seqs, reverseComplement(seqs))
  
  # --- 4. Sliding-window k-mer counts, collapsed across all sequences ---
  all_kmers <- mkAllStrings(c("A", "C", "G", "T"), k)
  results <- setNames(integer(length(all_kmers)), all_kmers)
  
  totals <- oligonucleotideFrequency(seqs2, width = k, step = 1,
                                     simplify.as = "collapsed")
  results[names(totals)] <- results[names(totals)] + totals
  
  # --- 5. Sort ascending and sum the middle 50% ---
  results <- sort(results, decreasing = FALSE)
  n <- length(results)
  middle <- results[(floor(n * 0.25) + 1):floor(n * 0.75)]
  middle_sum <- sum(middle)
  
  list(
    reads_orig            = n_orig,
    reads_after_length_filter = n_after_length_filtering,
    reads_after_N_filter  = n_after_ambiguity_filtering,
    middle_sum            = middle_sum
  )
}


lambda = function(sig_file, bg_file, expected_length, k = 8) {
  
  files <- c("sig" = sig_file,
             "bg"  = bg_file)
  df=do.call(rbind, lapply(files, function(f) 
    as.data.frame(count_kmer_middle(f, expected_length = expected_length, k=k))))
  
  N_sig=df["sig", "reads_after_N_filter"]
  N_bg=df["bg", "reads_after_N_filter"]
  B=df["bg", "middle_sum"]
  S=df["sig", "middle_sum"]
  
  lambda=(S/B) * (N_bg/N_sig)
  list(
    lambda = lambda,
    data=df,
    N_sig = N_sig,
    N_bg = N_bg,
    S = S,
    B = B
  )
}


if (sys.nframe() == 0) {
  option_list <- list(
    make_option(c("-s", "--sig"),    type = "character",
                help = "Signal file (.seq, .fasta, .fastq)"),
    make_option(c("-b", "--bg"),     type = "character",
                help = "Background file (.seq, .fasta, .fastq)"),
    make_option(c("-l", "--length"), type = "integer", default = 20,
                help = "Expected read length [default %default]"),
    make_option(c("-k", "--kmer"),   type = "integer", default = 8,
                help = "k-mer length [default %default]"),
    make_option(c("--header"), action = "store_true", default = TRUE,
                help = "Print TSV header line before the data row"),
    make_option(c("-o", "--out"),    type = "character", default = NULL,
                help = "Optional output TSV file for the data frame")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  if (is.null(opt$sig) || is.null(opt$bg)) {
    stop("Both --sig and --bg are required. Use --help for usage.")
  }
  
  result <- lambda(sig_file = opt$sig,
                   bg_file  = opt$bg,
                   expected_length = opt$length,
                   k = opt$kmer)
  
  
  
  # result <- lambda(sig_file = "/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/FOXO1_TTCAGC20NTA_AF_3.fastq.gz",
  #                  bg_file  = "/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/FOXO1_TTCAGC20NTA_AF_2.fastq.gz",
  #                  expected_length = 20,
  #                  k = 8)
  
  # Group      Total_sequences (After N filtering)  Hits Fraction Middle_8mer_total Max_8mer_count
  # 1 Background          151707                    272        0              1783634            171
  # 2 Signal              225026                    2440      1               2450095            269
  
  # $lambda
  # [1] 0.926
  
  #So the lambda would be 
  #(2450095* 151707)/( 1783634* 225026   )=0.9260834
  

  
  if (opt$header) {
    cat("sig_file\tbg_file\tN_sig\tN_bg\tS\tB\tlambda\n")
  }
  cat(sprintf("%s\t%s\t%d\t%d\t%d\t%d\t%.6f\n",
              basename(opt$sig), basename(opt$bg),
              result$N_sig, result$N_bg,
              result$S, result$B, result$lambda))
  
  if (!is.null(opt$out)) {
    write.table(result$data, opt$out, sep = "\t", quote = FALSE,
                col.names = NA)
    cat("Wrote", opt$out, "\n")
  }
}




