library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
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

library("readr")

file="/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_19032026.tsv" #3596

metadata=read_delim(file, delim="\t")

unique_combinations=unique(metadata[, c("CSC_SELEX_filename","CSC_SELEX_background_filename" )]) #2681
unique_combinations$ID=NA


for(i in 1:nrow(unique_combinations)){
  data_row=metadata %>% filter(CSC_SELEX_filename==unique_combinations$CSC_SELEX_filename[i] & CSC_SELEX_background_filename==unique_combinations$CSC_SELEX_background_filename[i])
  unique_combinations$ID[i]=data_row$ID[1]
}

unique_combinations=metadata %>% filter(ID %in% unique_combinations$ID)

for(i in 1:nrow(unique_combinations)){
  print(i)
  spacek=parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/",unique_combinations$ID[i],".txt"))
  
  unique_combinations$spacek_lambda[i]=spacek$lambda
  unique_combinations$spacek_N_sig[i]=spacek$summary %>% filter(Group=="Signal") %>% pull("Total_sequences")
  unique_combinations$spacek_N_bg[i]=spacek$summary %>% filter(Group=="Background") %>% pull("Total_sequences")
  unique_combinations$spacek_S[i]=spacek$summary %>% filter(Group=="Signal") %>% pull("Lower_half_8mer_total")
  unique_combinations$spacek_B[i]=spacek$summary %>% filter(Group=="Background") %>% pull("Lower_half_8mer_total")
    
}


#Read lambda etc. computed by the R script 
files=dir("/scratch/project_2013895/SELEX/lambda_computation/output") #2672

all_data=data.frame(signal_background_pair=gsub(".tsv", "", files), sig_file="", 
           bg_file="",  
           N_sig="",
           N_bg=NA,
           S=NA,
           B=NA,
           lambda=NA )

for(i in 1:length(files)){
  print(i)
  file=files[i]
  data=read_delim(paste0("/scratch/project_2013895/SELEX/lambda_computation/output/", file), delim="\t")[1,]

  all_data$sig_file[i]=data$sig_file
  all_data$bg_file[i]=data$bg_file
  all_data$N_sig[i]=data$N_sig
  all_data$N_bg[i]=data$N_bg
  all_data$S[i]=data$S
  all_data$B[i]=data$B
  all_data$lambda[i]=data$lambda
  
}


unique_combinations$signal_base=basename(unique_combinations$CSC_SELEX_filename)
unique_combinations$background_base=basename(unique_combinations$CSC_SELEX_background_filename)
unique_combinations$signal_background_pair=paste0(gsub(".fastq.gz", "",unique_combinations$signal_base), "_",
                                                  gsub(".fastq.gz", "",unique_combinations$background_base) )


unique_combinations=unique_combinations %>% left_join(all_data, by=c("signal_background_pair"))

#Filter out the corrupted files

unique_combinations$fastq_ftp_signal %>% is.na() %>% table()
unique_combinations$fastq_ftp_background %>% is.na() %>% table()

unique_combinations=unique_combinations %>% filter(!is.na(fastq_ftp_signal)) #2672

setwd("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments")

#saveRDS(unique_combinations,"RData/lambdas.Rds")
unique_combinations = readRDS("RData/lambdas.Rds")

abs(unique_combinations$spacek_lambda-unique_combinations$lambda) %>% max(na.rm=TRUE) #5e-04 0.0005
abs(unique_combinations$spacek_lambda-unique_combinations$lambda) %>% min(na.rm=TRUE) #0


(unique_combinations$spacek_N_sig==unique_combinations$N_sig) %>% table() #All true
(unique_combinations$spacek_N_bg==unique_combinations$N_bg) %>% table() #All true
(unique_combinations$spacek_S==unique_combinations$S) %>% table() # all true
(unique_combinations$spacek_B==unique_combinations$B) %>% table() #All true

library(scales)

ggplot(unique_combinations, aes(x = lambda)) +
  geom_histogram(bins = 60,
                 fill = "steelblue",   # inside color
                 color = "black" ) +   # border color
  scale_x_continuous(labels = comma, 
                     breaks = seq(0, 1.2, by = 0.2)) +
  labs(
    x = "Lambda",
    y = "Count",
    title = "Distribution of Lambdas"
  ) +
  theme_classic()


unique_combinations %>%filter(lambda>1) %>% nrow() #26

