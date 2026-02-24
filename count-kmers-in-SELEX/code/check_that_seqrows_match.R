library("readr")
library("dplyr")
library("tidyverse")
library(R.utils)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(readr)
library(dplyr)

metadata=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv", delim="\t")


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



metadata$signal_fastq_rownro=NA
metadata$signal_fasta_rownro=NA
metadata$signal_seq_rownro=NA
metadata$signal_fasta_seqlength=NA

metadata$background_fastq_rownro=NA 
metadata$background_fasta_rownro=NA  
metadata$background_seq_rownro=NA
metadata$background_fasta_seqlength=NA 

metadata$spacek_signal_sequences=NA
metadata$spacek_background_sequences=NA
metadata$lambda=NA

metadata$kmer_counts_exist=NA

for( i in 1:nrow(metadata) ){
  print(i)
  metadata$signal_fastq_rownro[i] <- as.integer(system(paste0("zcat ",metadata$CSC_SELEX_filename[i]," | wc -l"), intern = TRUE))
  metadata$background_fastq_rownro[i] <- as.integer(system(paste0("zcat ",metadata$CSC_SELEX_background_filename[i]," | wc -l"), intern = TRUE))
  
  seqkit_sample <- try(read_table(
    paste0("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/seqkit_stats/",
           sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_filename[i])), ".txt"),
    col_types = cols(
      file     = col_character(),
      format   = col_character(),
      type     = col_character(),
      num_seqs = col_character(),
      sum_len  = col_character(),
      min_len  = col_integer(),
      avg_len  = col_double(),
      max_len  = col_integer()
    )
  ), silent=TRUE)



  if(!(inherits(seqkit_sample,'try-error'))& nrow(seqkit_sample)!=0){
    metadata$signal_fasta_rownro[i] <- as.integer(gsub(",","",seqkit_sample$num_seqs))
    metadata$signal_fasta_seqlength[i]=unique(c(seqkit_sample$min_len, seqkit_sample$avg_len, seqkit_sample$max_len)[!is.na(c(seqkit_sample$min_len, seqkit_sample$avg_len, seqkit_sample$max_len))])
  }
  seqkit_background <- try(read_table(
    paste0("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/seqkit_stats/",
           sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_background_filename[i])), ".txt"),
    col_types = cols(
      file     = col_character(),
      format   = col_character(),
      type     = col_character(),
      num_seqs = col_character(),
      sum_len  = col_character(),
      min_len  = col_integer(),
      avg_len  = col_double(),
      max_len  = col_integer()
    )
  ), silent=TRUE)

  if(!(inherits(seqkit_background,'try-error')) & nrow(seqkit_background)!=0){
    metadata$background_fasta_rownro[i] <- as.integer(gsub(",","",seqkit_background$num_seqs))
    metadata$background_fasta_seqlength[i]=unique(c(seqkit_background$min_len, seqkit_background$avg_len, seqkit_background$max_len)[!is.na(c(seqkit_background$min_len,
                                                                                                                                    seqkit_background$avg_len,
                                                                                                                                           seqkit_background$max_len))])
  }

  if(file.exists(paste0("/scratch/project_2013895/SELEX/preprocessed_data/",
                        sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_filename[i])), ".seq"))){
    metadata$signal_seq_rownro[i]=signal_seq_rownro=try( as.integer(sub("^\\s*(\\d+).*", "\\1",
                                                                        system(paste0("wc -l ", paste0("/scratch/project_2013895/SELEX/preprocessed_data/",
                                                                                                       sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_filename[i])), ".seq")
                                                                        ), intern = TRUE))), silent=TRUE)
  }



  if(file.exists(paste0("/scratch/project_2013895/SELEX/preprocessed_data/",
                        sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_background_filename[i])), ".seq"))){
  metadata$background_seq_rownro[i]=as.integer(sub("^\\s*(\\d+).*", "\\1",
                                               system(paste0("wc -l ", paste0("/scratch/project_2013895/SELEX/preprocessed_data/",
                                                                              sub("\\.(fastq|fq)(\\.gz)?$", "", basename(metadata$CSC_SELEX_background_filename[i])), ".seq")
                                               ), intern = TRUE)))

  }

  metadata$kmer_counts_exist[i]=file.exists(paste0("/scratch/project_2013895/SELEX/spacek/kmer_counts/", metadata$ID[i], ".txt"))

  #Read signal and background from spacek output and lambda also
  ## ---- usage ----
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", metadata$ID[i], ".txt")), silent=TRUE)
  if(class(res)=="try-error"){
    next
  }
  #res$cutoff      # numeric
  #res$files       # data.frame with file_id, eof_line, num_sequences
  #res$summary     # tibble with Background / Signal rows
  #res$lambda      # numeric

  metadata$spacek_signal_sequences[i]=res$summary %>% filter(Group=="Signal") %>% pull(Total_sequences)
  metadata$spacek_background_sequences[i]=res$summary %>% filter(Group=="Background") %>% pull(Total_sequences)
  metadata$lambda[i]=res$lambda


  
}

write_delim(metadata, "/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_with_seqinfo_spacek_output.tsv", delim="\t")



# Analysis of the results -------------------------------------------------

metadata=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_with_seqinfo_spacek_output.tsv", delim="\t")

metadata$signal_fastq_rownro[1]
metadata$signal_fasta_rownro[1]
metadata$signal_fasta_rownro[1]*4
metadata$signal_seq_rownro[1]
metadata$signal_seq_rownro[1]*4
metadata$spacek_signal_sequences[1]*4

(metadata$signal_seq_rownro*4==metadata$signal_fastq_rownro) %>% table()
#FALSE  TRUE 
#8  3267 

test=metadata[which(!(metadata$signal_seq_rownro*4==metadata$signal_fastq_rownro)),]

inds=which(!(metadata$signal_seq_rownro*4==metadata$signal_fastq_rownro))
inds[1]
metadata$signal_fastq_rownro[inds[2]]
metadata$signal_fasta_rownro[inds[2]]
metadata$signal_fasta_rownro[inds[2]]*4
metadata$signal_seq_rownro[inds[2]]
metadata$signal_seq_rownro[inds[2]]*4
metadata$spacek_signal_sequences[inds[2]]*4

(metadata$background_seq_rownro*4==metadata$background_fastq_rownro) %>% table()
#FALSE  TRUE 
#61  3366 

(metadata$signal_seq_rownro==metadata$spacek_signal_sequences) %>% table()
#FALSE  TRUE 
#2963    19 

(metadata$background_seq_rownro==metadata$spacek_background_sequences) %>% table()
#FALSE  TRUE 
#2946    36 

(!is.na(metadata$lambda)) %>% table()
#FALSE  TRUE 
#635  2938 

metadata$kmer_counts_exist %>% table(useNA = "always")

metadata$signal_seq_rownro %>% hist()

#Do we have the k-mers with motif match scores and SELEX scores 

results_path="/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/"
files=dir(results_path)

motif_match_files=files[-grep("_SELEX_counts.tsv",files)]
motif_match_files=gsub(".tsv", "", motif_match_files)

SELEX_files=files[grep("_SELEX_counts.tsv",files)]
SELEX_files=gsub("_SELEX_counts.tsv", "", SELEX_files)


metadata$motif_match_scores=metadata$ID %in% motif_match_files
metadata$SELEX_scores=metadata$ID %in% SELEX_files

FOXI1_HT-SELEX_TCGGAA20NGA_AF_GTAAACA_1_4

#metadata$ID[i]
file <- paste0("/scratch/project_2013895/SELEX/spacek/kmer_counts/", "MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6", ".txt")

# detect header line
all_lines <- read_lines(file)
hdr_idx <- which(startsWith(all_lines, "Background\tSignal\tSequence_length\tGap_position"))[1]
cols <- strsplit(all_lines[hdr_idx], "\t")[[1]]

# read the data rows after the header; keep exactly header’s columns
df <- read_tsv(file, skip = hdr_idx, col_names = FALSE, show_col_types = FALSE)
names(df) <- c(cols[1:(length(cols)-1)], "local_max_type", "value1", cols[length(cols)], "value2")

#remove rows from the end 

df=df[-which(df$Background!=df$Background[1]),]

df2 = df %>%
  separate(
    col = value2,
    into = c("value2", "local_hamming"),
    sep = "(?:\\t|\\\\t|\\s+)",  # tab, literal "\t", or whitespace
    extra = "merge",             # keep anything after the first split together
    fill  = "right"              # if missing flag, put NA
  ) %>%
  mutate(value2 = readr::parse_number(value2))


test=df2 %>% filter(!is.na(local_hamming))

#df2 has 16384 rows, seed length is 7, the number of possible k-mers is 4^7 16384

df2$Kmer %>% unique() %>% length() #16384

#Load the seed realisations and their Hamming 1 distance k-mers
kmers_Hamming1=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",
                                 "MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6", 
                                 ".tsv"), col_names=FALSE) #224
names(kmers_Hamming1)=c("Kmer", "score", "rank")

#Load the seed realisations and their Hamming 1&2 distance k-mers
kmers_Hamming12=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",
                                 "MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6", "_score_rank", 
                                 ".tsv"), col_names=FALSE) #1280
names(kmers_Hamming12)=c("Kmer", "score", "rank")
#Load all canonical k-mers of length 7

kmers_all=read_delim("/scratch/project_2006203/TFBS/Results/kmers/kmers_k7.canonical_cpp.txt", delim=" ",col_names=FALSE)
names(kmers_all)="Kmer" #8192, this is 16384/2

#find . -maxdepth 1 -type f -printf '%s\t%f\n' | sort -n | numfmt --to=iec --field=1 --padding=7 | column -t


match_ind=match(kmers_Hamming12$Kmer, df2$Kmer)

head(kmers_Hamming12$Kmer)
head(df2$Kmer[match_ind])

test=left_join(kmers_Hamming12, df2 %>% select("Kmer", "Signal count", "Background count", "Fold_change", "Local_max", "local_max_type", "local_hamming"), by = "Kmer")

test2=left_join(kmers_all, df2 %>% select("Kmer", "Signal count", "Background count", "Fold_change", "Local_max", "local_max_type", "local_hamming"), by = "Kmer")
