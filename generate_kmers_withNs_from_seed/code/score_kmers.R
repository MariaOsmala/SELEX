library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("tidyverse")
library(R.utils)
library(stringr)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)
arrays=as.numeric(args[1])
addition=as.numeric(args[2])
arrays=arrays+addition

#length=10 
length=1 #starting from 231 2301-2310
start_ind=(arrays-1)*length+1 #100
end_ind=arrays*length #100


metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635

#which(metadata$seed=="NNNACGANNNNNNTCGTNNN") #992
#which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3")
if(end_ind > nrow(metadata)){
  end_ind=nrow(metadata)
}

#sort based on the seed length

metadata=metadata[order(metadata$length),]

#which(metadata$ID=="Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3") #373
#which(metadata$ID=="BARHL2_HT-SELEX_TCCAGT40NGAC_AI_NNTAATTGNN_1_3") #377
#which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3") #2333
results_path="/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/"

#results=dir(results_path)
#results=gsub(".tsv", "", results) #3461
#missing=metadata %>% filter(!(ID %in% results)) %>% pull(ID) #175
#missing=metadata %>% filter(!(ID %in% results))
#missing %in% failed_experiments_all$motif %>% table()
#FALSE  TRUE 
#3484   175 
#missing %in% failed_experiments$motif %>% table()
#FALSE  TRUE 
#29   146 
#missing=missing[!(missing %in% failed_experiments$motif)]
#there are some additional failed cases
#failed_experiments_all %>% filter(motif %in%missing) %>% View()

path_local_scratch=Sys.getenv("LOCAL_SCRATCH")
#dir(path_local_scratch)

lambda_path="/scratch/project_2013895/SELEX/spacek/lambda/"
# ls | wc -l 3573

#For which motifs lambda fails?

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

  
#index=arrays

for(index in seq(start_ind, end_ind, 1)){ #0-9

  
  print(index)
  data=metadata[index,] 
  seed=data$seed
    
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", data$ID, ".txt")), silent=TRUE)
  
  if(class(res)!="try-error"){
    N_sig=res$summary %>% filter(Group=="Signal") %>% pull(Total_sequences)
    N_bg=res$summary %>% filter(Group=="Background") %>% pull(Total_sequences)
    lambda=res$lambda
  
  # Score k-mers containing Ns based on k-mer counts 
  
  
  kmer_signal_count_file=paste0("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/", data$ID, "_signal.txt")
  kmer_background_count_file=paste0("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/", data$ID, "_background.txt")
  
  kmers_signal_counts=read_delim(kmer_signal_count_file, col_names=FALSE, delim="\t") #1126
  kmers_background_counts=read_delim(kmer_background_count_file, col_names=FALSE, delim="\t") #1126
  
  names(kmers_signal_counts)=c("Kmer","signal_count")
  names(kmers_background_counts)=c("Kmer","background_count")
  
  counts=left_join(kmers_signal_counts, kmers_background_counts, by="Kmer")
  
               
  counts$`Corrected_count`=counts$signal_count-(lambda*(N_sig/N_bg))*counts$background_count
  
  #Order based on the Corrected count, divide by maximum count, and add rank
  
  counts <- counts %>%
    mutate(`Corrected_count_scaled`=`Corrected_count`/max(`Corrected_count`)) %>% 
    arrange(desc(`Corrected_count_scaled`)) %>%
    mutate(`Corrected_count_rank` = row_number()) %>% select(Kmer, `Corrected_count_scaled`, Corrected_count_rank)
  
  write_delim(counts, #244
              file=paste0(results_path,data$ID,".tsv"),
              delim="\t")
    

  }else{
    print("no lambda")
  } #lambda exists
}

  