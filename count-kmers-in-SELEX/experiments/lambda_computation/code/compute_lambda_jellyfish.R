#!/usr/bin/env Rscript


# Functions ---------------------------------------------------------------


# --- Count reads (FASTQ: 4 lines per record) ---
count_reads <- function(fq) {
  if (grepl("\\.gz$", fq)) {
    n <- as.integer(system2("zcat", c(fq, "| wc -l"), stdout = TRUE))
  } else {
    n <- as.integer(system2("wc", c("-l", fq), stdout = TRUE))
  }
  n / 4
}


# --- Load counts computed by jellyfish and compute IQR sum ---
library(Biostrings)
library("dplyr")
load_counts <- function(path) {
  df <- read.table(path, header = FALSE, col.names = c("kmer", "count")) #65536
  
  k <- 8
  bases <- c("A", "C", "G", "T")
  all_kmers <- do.call(paste0, expand.grid(rep(list(bases), k)))
  #only unique
  df2=data.frame(kmer=unique(all_kmers)) #65536
  
  df2=df2 %>% left_join(df, by="kmer")
  df2$count[is.na(df2$count)]=0
  sort(df2$count)
}

iqr_sum <- function(counts) {
  n  <- length(counts) #32886
  lo <- floor(n / 4) + 1
  hi <- floor(3 * n / 4)
  sum(counts[lo:hi])
}

# --- Parse spacek reports ---
library(stringr)
parse_spacek_report <- function(path) {
  x <- read_lines(path)
  
  ## cutoff
  cutoff <- x[str_detect(x, "^cutoff:")]
  cutoff <- as.numeric(str_match(cutoff, "^cutoff:\\s*([+-]?[0-9.]+)")[,2])
  
  ## EOF encountered + Number of sequences (pair per file)
  eof_idx  <- which(str_detect(x, "^EOF encountered in file\\s+\\d+\\s+on line\\s+\\d+"))
  num_idx  <- which(str_detect(x, "^Number of sequences in file\\s+\\d+:\\s+\\d+"))
  
  eof_dat <- tibble(line = x[eof_idx]) |>
    mutate(file_id = as.integer(str_match(line, "file\\s+(\\d+)")[,2]),
           eof_line = as.integer(str_match(line, "line\\s+(\\d+)")[,2])) |>
    select(file_id, eof_line)
  
  num_dat <- tibble(line = x[num_idx]) |>
    mutate(file_id = as.integer(str_match(line, "file\\s+(\\d+)")[,2]),
           num_sequences = as.integer(str_match(line, ":(\\s*\\d+)")[,2])) |>
    select(file_id, num_sequences)
  
  files <- eof_dat |> full_join(num_dat, by = "file_id") |> arrange(file_id)
  
  ## Small 2-row summary table (Background / Signal)
  # find header line that starts with a tab then "Total sequences"
  hdr_i <- which(str_detect(x, "^\\s*Total sequences\\s+Hits\\s+Fraction\\s+Lower half 8mer total\\s+Max 8mer count"))[1]
  if (is.na(hdr_i)) stop("Couldn't find the Total sequences table header.")
  
  # Next two non-empty lines are Background and Signal (strip possible extra spaces)
  table_lines <- x[(hdr_i+1):(hdr_i+3)]
  table_lines <- table_lines[nzchar(str_trim(table_lines))]  # drop blank
  # Build a mini text block: header + rows
  mini <- c(
    "Group\tTotal_sequences\tHits\tFraction\tLower_half_8mer_total\tMax_8mer_count",
    str_replace_all(str_trim(table_lines), "\\s+", "\t")
  )
  summary_tbl <- read_tsv(I(mini), show_col_types = FALSE,
                          col_types = cols(
                            Group = col_character(),
                            Total_sequences = col_integer(),
                            Hits = col_integer(),
                            Fraction = col_double(),
                            Lower_half_8mer_total = col_integer(),
                            Max_8mer_count = col_integer()
                          ))
  
  ## Lambda
  lambda_line <- x[str_detect(x, "^Lambda\\s+")]
  lambda <- as.numeric(str_match(lambda_line, "^Lambda\\s+([0-9.]+)")[,2])
  
  list(
    cutoff = cutoff,
    files  = files,
    summary = summary_tbl,
    lambda = lambda
  )
}

# Main --------------------------------------------------------------------
library("readr")

file="/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_19032026.tsv" #3596

metadata=read_delim(file, delim="\t")

#grep("FOXO1_TTCAGC20NTA_AF_3", metadata$CSC_SELEX_filename)
i=423

#Before N filtering
#225057
N_sig <- count_reads(metadata$CSC_SELEX_filename[i])
#151717                     
N_bg  <- count_reads(metadata$CSC_SELEX_background_filename[i])


#After N filtering
#225026
N_sig <- count_reads(paste0("/scratch/project_2013895/SELEX/data_exclude_reads_with_N/", 
                            basename(metadata$CSC_SELEX_filename[i]))
                          )
#151707                     
N_bg  <- count_reads(paste0("/scratch/project_2013895/SELEX/data_exclude_reads_with_N/", 
                            basename(metadata$CSC_SELEX_background_filename[i]))
)




sig_counts <- load_counts(paste0("/scratch/project_2013895/SELEX/lambda_computation/filtered_counts_without_C/", 
                                 gsub(".fastq.gz","",basename(metadata$CSC_SELEX_filename[i])), ".txt")
)
bg_counts  <- load_counts(paste0("/scratch/project_2013895/SELEX/lambda_computation/filtered_counts_without_C/", 
                                 gsub(".fastq.gz","",basename(metadata$CSC_SELEX_background_filename[i])), ".txt")
)

S <- iqr_sum(sig_counts) #2450095
B <- iqr_sum(bg_counts) #1783634

lambda <- (S * N_bg) / (B * N_sig) #0.9260834

message(sprintf("N_sig: %d", N_sig)) #225026
message(sprintf("N_bg:  %d", N_bg)) #151707
message(sprintf("S (signal IQR sum):     %d", S))
message(sprintf("B (background IQR sum): %d", B))
message(sprintf("lambda = %.6f", lambda))

parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/",metadata$ID[i],".txt"))

# Group      Total_sequences (After N filtering)  Hits Fraction Middle_8mer_total Max_8mer_count
# 1 Background          151707                    272        0              1783634            171
# 2 Signal              225026                    2440      1               2450095            269

# $lambda
# [1] 0.926

#So the lambda would be 
#(2450095* 151707)/( 1783634* 225026   )=0.9260834



# Spacek output
# cutoff: -0.000100: This is the PWM match threshold used in your run. 
# In multinomial mode it is set as -m - 0.0001; with m=0 it becomes -0.0001 at spacek40.c:9703.

# EOF encountered in file 0 on line 107035
# 
# Number of sequences in file 0: 107034
# EOF encountered in file 1 on line 41580
# 
# Number of sequences in file 1: 41579

#Column labels like -9 ... 29 are positions relative to the seed anchor window.
# All Hits: counts of bases at each aligned position for all detected hits in the signal set.
# All Hits        Position
# All Hits                -9      -8      -7      -6      -5      -4      -3      -2      -1      0       1       2       3       4       5       6       7       8   910      11      12      13      14      15      16      17      18      19      20      21      22      23      24      25      26      27      28      29
# All Hits        A       0       0       0       0       0       0       0       0       0       62      0       0       0       0       0       0       0       0   106      0       0       110     110     0       0       43      70      67      0       0       0       0       0       0       0       0       0       0       0
# All Hits        C       0       0       0       0       0       0       0       0       0       53      0       0       141     0       0       0       146     0   00       0       0       0       0       0       92      55      0       0       0       0       0       0       0       0       0       0       0       0
# All Hits        G       0       0       0       0       0       0       0       0       25      0       0       71      0       0       0       0       0       0   066      65      0       0       0       0       0       0       26      43      0       0       0       0       0       0       0       0       0       0
# All Hits        T       0       0       0       0       0       0       0       0       0       0       80      0       0       81      81      82      0       81  00       0       0       0       82      80      0       0       0       0       0       0       0       0       0       0       0       0       0       0
#
# counts from sequences that have exactly one hit.
# One Hit         Position
# One Hit                 -9      -8      -7      -6      -5      -4      -3      -2      -1      0       1       2       3       4       5       6       7       8   910      11      12      13      14      15      16      17      18      19      20      21      22      23      24      25      26      27      28      29
# One Hit         A       0       0       0       0       0       0       0       0       0       62      0       0       0       0       0       0       0       0   106      0       0       110     110     0       0       43      70      67      0       0       0       0       0       0       0       0       0       0       0
# One Hit         C       0       0       0       0       0       0       0       0       0       53      0       0       141     0       0       0       146     0   00       0       0       0       0       0       92      55      0       0       0       0       0       0       0       0       0       0       0       0
# One Hit         G       0       0       0       0       0       0       0       0       25      0       0       71      0       0       0       0       0       0   066      65      0       0       0       0       0       0       26      43      0       0       0       0       0       0       0       0       0       0
# One Hit         T       0       0       0       0       0       0       0       0       0       0       80      0       0       81      81      82      0       81  00       0       0       0       82      80      0       0       0       0       0       0       0       0       0       0       0       0       0       0
#
# This table is only populated in a specific branch requiring multinomial background correction conditions at spacek40.c:10354.
# OneHit_Exp      Position
# OneHit_Exp              -9      -8      -7      -6      -5      -4      -3      -2      -1      0       1       2       3       4       5       6       7       8   9   10      11      12      13      14      15      16      17      18      19      20      21      22      23      24      25      26      27      28      29
# OneHit_Exp      A       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   0   0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
# OneHit_Exp      C       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   0   0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
# OneHit_Exp      G       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   0   0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
# OneHit_Exp      T       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   0   0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
# Total number of two hits: 0
# Total possible spacings: 0
# 
# ! Two or more hits connecting 
# Two-hit statistics
# 
# Total number of two hits: 0 and Total possible spacings: 0
# No exact two-hit events were collected in the dedicated two-hit PWM accumulator spacek40.c:10572.
# Two or more hits connecting matrix has 92 at Head-to-tail spacing 0:
# This looks surprising but comes from how that matrix is filled: it includes self-pairs (counter2 starts from counter), so each hit contributes a spacing-0 head-to-tail count spacek40.c:10139. That is why you see 92 there even though exact two-hit matrix is all zeros.
# Exactly two hits connecting matrix all zeros:
# This is the strict two-hit-only matrix, and it confirms no true two-hit configurations were observed spacek40.c:10642.

# ! Orientation                   Spacing
# !                       0       1       2       3
# ! ------------------------------------------------------------------------------------------------------------------------
# ! Head to tail > > :    92      0       0       0
# ! Head to head > < :    0       0       0       0
# ! Tail to tail < > :    0       0       0       0
# 
# 2! Exactly two hits connecting matrix
# 2! Orientation                  Spacing
# 2!                      0       1       2       3
# 2! ------------------------------------------------------------------------------------------------------------------------
# 2! Head to tail > > :   0       0       0       0
# 2! Head to head > < :   0       0       0       0
# 2! Tail to tail < > :   0       0       0       0
# 
# +
# Hits are strongly position-biased to reverse-strand positions 2 and 3 in the allowed window.
# Information content 2.32 (forward) and 1.04 (reverse) quantifies positional concentration (higher means more concentrated/non-uniform), computed at spacek40.c:10938.
# + Hit positions 1       2       3       4       Information content
# + Forward:      0       0       0       0       2.32
# + Reverse:      0       35      57      0       1.04
#
# Background vs signal summary
# Background: 107034 sequences, 3 hit-containing
# Signal: 41579 sequences, 92 hit-containing
# Fraction prints as 0 for both because formatting is %.0f (rounded whole percent) at spacek40.c:10950.
# Actual fractions are about 0.003% (background) and 0.221% (signal), so signal is clearly enriched.
# Lambda 0.871
# Size normalization factor computed from middle (25-75%) kmer count mass and sequence totals spacek40.c:10332.

# the exact print and calculation points in spacek40.c:10313 and spacek40.c:10950
# The program takes all 8-mer counts (or k-mer length set by -kl), sorts them, 
# and sums the middle 50% of counts (from 25th to 75th percentile), not literally the lowest half.
# It does this separately for background and signal, then uses those sums to compute lambda.
# Code path: spacek40.c:10313, spacek40.c:10332, and printed at spacek40.c:10950.

# Total sequences Hits    Fraction        Lower half 8mer total   Max 8mer count
# Background      107034 (N_bg)  3       0       1213416(B) 233
# Signal          41579 (N_sg)   92      0       410443(S)  114
# Lambda 0.871
# 
#So the lambda would be 
#(410443*107034)/(1213416*41579)=0.871

#Seed output
#            Total sequences	Hits	Fraction	Lower half 8mer total	Max 8mer count
#Background	63037 (N_bg)	          1921	3	        742948 (B)	              424
#Signal    	80122	 (N_sg)          11039	14	      838913 (S)	              3292
#Lambda 0.888

#So the lambda would be 
#(838913*63037)/(742948*80122)=0.8883871


# the PWM/profile is most similar to the internal reference 6-mer associated with HSF in the built-in TF-kmer catalog, using the scoring routine around spacek40.c:2423 and print at spacek40.c:2541.
# Best characteristic kmer match for All hits background  : 3 is number 51: HSF with score 6.00 at position 29 strand 1
# 
# Best characteristic kmer match for All hits uncorrected  : 92 is number 51: HSF with score 6.00 at position 29 strand 1
# 
# Best characteristic kmer match for All hits : 92 is number 51: HSF with score 6.00 at position 29 strand 1
# 
# Best characteristic kmer match for One hit   : 92 is number 51: HSF with score 6.00 at position 29 strand 1
#
# Because that PWM was all zeros, matcher stays at default worst score baseline (-100), so it reports NONE.
# Best characteristic kmer match for One hit exponential : 92 is number 0: NONE with score -100.00 at position 0 strand 0
# 
# Best characteristic kmer match for Forward backgr is number 59: ZNF143 with score 1.98 at position 8 strand 1
# 
# Best characteristic kmer match for Reverse backgr is number 59: ZNF143 with score 1.99 at position 19 strand -1
# 
# Signal 0.66 bits vs background 0.19 bits for all sequences: signal kmer spectrum is more structured/less uniform.
# “Sequences not hit” values are almost the same as “All sequences”, meaning removing hit-containing sequences barely changes 
# global 8-mer complexity (expected since hit-containing sequences are a small fraction).
# Information content for kmer length 8   All sequences   Sequences not hit
# Signal                                  0.66    0.65    bits
# Background                              0.19    0.19    bits
#
# no palindromic forward/reverse duplicate-hit cases were removed under the active palindrome handling logic.
# Palindromic hits
# Background:     0
# Signal    :     0
# 
# 
# Type    %_of_all        %_of_1or2       OrSeed  AndSeed PWMSeed
# OneHit  100.0%  100.0%  TGCTTTCTAGGAATTMM       TGCTTTCTAGGAATTMM       MTGCTTTCTAGGAATT        onehit_seed
# 
# ** user specified seed TGCTTTCTAGGAATTMM with no seed iteration
# 
# ** PWM SEED CALCULATION
# * full seed from PWM: ??????????????????GMTGCTTTCTAGGAATTMMRG???
#   * trimmed to round max length+2: GMTGCTTTCTAGGAATTM
# * end Ns removed: GMTGCTTTCTAGGAATTM
# * seed GMTGCTTTCTAGGAATTM is too long, removes last base M
# * seed GMTGCTTTCTAGGAATT is too long, removes first base G
# 
# version spacek40 v0.200 DEC 19 2018 : command spacek40 -20N _-f 
# /scratch/project_2013895/SELEX/preprocessed_data/BCL6B_TGCGGG20NGA_AC_3.seq /scratch/project_2013895/SELEX/preprocessed_data/BCL6B_TGCGGG20NGA_AC_4.seq TGCTTTCTAGGAATTMM 


