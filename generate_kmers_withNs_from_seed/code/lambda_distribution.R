library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("tidyverse")
library(R.utils)
library(stringr)
library(tidyr)

rm(list=ls())


metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX")) #3636
metadata=metadata %>% filter(!is.na(seed)) #3635

metadata=metadata[order(metadata$length),]

results_path="/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/Data/"

lambda_path="/scratch/project_2013895/SELEX/spacek/lambda/"
#Remove first zero-sized files
#find . -maxdepth 1 -type f -size 0 -delete
# ls | wc -l 3470

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

metadata$lambda=NA
metadata$N_sig=NA
metadata$N_bg=NA


for(index in 1:nrow(metadata)){ #0-9

  #index=1
  print(index)
  data=metadata[index,] 
  #data=tmp[2,]
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", data$ID, ".txt")), silent=TRUE)
  
  if(class(res)!="try-error"){
    #lambda exists
    N_sig=res$summary %>% filter(Group=="Signal") %>% pull(Total_sequences)
    N_bg=res$summary %>% filter(Group=="Background") %>% pull(Total_sequences)
    lambda=res$lambda
    metadata$lambda[index]=lambda
    metadata$N_sig[index]=N_sig
    metadata$N_bg[index]=N_bg
    
   }else{
    print("no lambda")
    
  }
}


metadata$lambda %>% is.na() %>% table()
#FALSE  TRUE 
#FALSE  TRUE 
#3596    39 

tmp=metadata %>% filter(is.na(lambda)) #%>% pull(ID) %>% sort() %>% head()

#For how many motifs, either signal or background is missing
table(is.na(metadata$CSC_SELEX_filename) | is.na(metadata$CSC_SELEX_background_filename)) #39

#for how many motifs, signal and background exist but lambda computation fails, none
table( (!is.na(metadata$CSC_SELEX_filename) & !is.na(metadata$CSC_SELEX_background_filename)) & is.na(metadata$lambda) )

ggplot(metadata, aes(x = lambda)) +
  geom_histogram(#aes(y = after_stat(density)),
    na.rm = TRUE,
                 binwidth=0.01,fill = "steelblue",
                 color = "black")+
  #geom_density(color = "red", linewidth = 1)+
  theme_minimal()+
  coord_cartesian(xlim = c(0, 1.2)) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1))
  
  

table(metadata$lambda >1) #25

#Indicate for which motifs we got scored k-mers

kmer_scores=dir("/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds_June2026/") #3596

kmer_scores=gsub(".tsv", "", kmer_scores)

metadata$kmer_scoring_exists=FALSE
metadata$kmer_scoring_exists[which(metadata$ID %in% kmer_scores)]=TRUE

metadata %>% select(kmer_scoring_exists) %>% table() #FALSE 40

write_delim(metadata,
            file="/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026_lambda_info.tsv",
            delim="\t")


#metadata=read_delim(file="/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026_lambda_info.tsv",
#            delim="\t")
  