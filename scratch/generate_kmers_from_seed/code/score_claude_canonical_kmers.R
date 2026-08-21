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

if(end_ind > nrow(metadata)){
  end_ind=nrow(metadata)
}

#sort based on the seed length

metadata=metadata[order(metadata$length),]

#which(metadata$ID=="Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3") #373
#which(metadata$ID=="BARHL2_HT-SELEX_TCCAGT40NGAC_AI_NNTAATTGNN_1_3") #377

results_path="/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/"

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

pseudocount=0.01

path_local_scratch=Sys.getenv("LOCAL_SCRATCH")
#dir(path_local_scratch)

# Function to create PWM from PFM -----------------------------------------

#C is PFM, p is pseudocount
PPMp <- function(C, p){
  (C + p / nrow(C)) / matrix( (colSums(C) + p),nrow=4, ncol=ncol(C),byrow=TRUE)
  
}

# Position-specific scoring matrix ----------------------------------------

#Position specific scoring matrix
PSSM <- function(C, B) log(C / B) #ln

score_sequence <- function(seq, pssm){
  #seq="CAGTGTGGTCGC"
  seq=sapply(as.character(seq), strsplit, "")[[1]]
  sum(as.numeric(mapply(function(base, pos) pssm[base, pos], seq, seq_along(seq))))
}

score_kmers_with_motif <- function(kmer_file, pssm, output_file) {
  in_con <- file(kmer_file, open = "r")
  out_con <- file(output_file, open = "w")
  on.exit({ close(in_con); close(out_con) })
  
  while (length(line <- readLines(in_con, n = 1, warn = FALSE)) > 0) {
    
    kmer_dna <- DNAString(line)
    rc_kmer_dna <- reverseComplement(DNAString(kmer_dna))
    
    # I need to score both the k-mer and its reverse complement and take the best score
    # Write kmer and score
    writeLines(paste(kmer_dna, 
                     max(score_sequence(as.character(kmer_dna), pssm),
                         score_sequence(as.character(rc_kmer_dna), pssm)
                     ), 
                     sep = "\t"), out_con)
  }
}
  

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
    
  # Score the k-mers by pwm -------------------------------------------------
  
  pcm=as.matrix(read.table(file=gsub("../../","/projappl/project_2006203/TFBS/",data$filename)))
  dimnames(pcm)=list(c("A", "C", "G", "T"))
    
  pfm <- TFBSTools::PFMatrix(ID=data$ID, 
                               strand="+", 
                               bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                               profileMatrix=pcm
    )
    
    
  #Check that two approaches gives the same result
  pwm=TFBSTools::toPWM(pfm, type="prob", pseudocounts=pseudocount, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  PWM=PPMp(pfm@profileMatrix,pseudocount)
    
  pssm=PSSM(PWM, 0.25)
  
  
  #real_file=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds","/",data$ID,"_canonical_realizations.txt")
  #hamming1_file=paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds","/",data$ID,"_canonical_hamming1.txt")
  
  all_kmers=paste0("/scratch/project_2013895/SELEX/combined_kmers_from_long_degenerate_seeds","/",data$ID,".txt")
  #test=read_table(all_kmers, col_names=FALSE) #224
  
  #file.create(all_kmers)
  
  #file.append(all_kmers, 
  #            c(real_file, hamming1_file))
  
  output_file=paste0(path_local_scratch,"/",data$ID,"_score.tsv")
    

  
  score_kmers_with_motif(kmer_file=all_kmers, pssm=pssm, output_file=output_file) 
  #dir(path_local_scratch)  
  #test=read_table(output_file, col_names=FALSE) #224
  
  output_file_score_rank=paste0(path_local_scratch,"/",data$ID,".tsv")
    
  # Construct the shell command
  cmd <- sprintf("sort -k2,2nr %s | awk '{printf \"%%s\\t%%s\\t%%d\\n\", $1, $2, NR}' > %s", output_file, output_file_score_rank)
    
  # Run it
  system(cmd)
  system(paste0("rm ",output_file))
  system(paste0("cp ", output_file_score_rank," ",results_path ))
  system(paste0("rm ",output_file_score_rank))
  
  #test=read_table(paste0(path_local_scratch,"/",data$ID,".tsv"), col_names=FALSE)  
  #Score based on SELEX data 
  
  #extract lambda if possible
  
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", data$ID, ".txt")), silent=TRUE)
  
  if(class(res)!="try-error"){
    N_sig=res$summary %>% filter(Group=="Signal") %>% pull(Total_sequences)
    N_bg=res$summary %>% filter(Group=="Background") %>% pull(Total_sequences)
    lambda=res$lambda
  
  
  kmer_count_file=paste0("/scratch/project_2013895/SELEX/spacek/kmer_counts/", data$ID, ".txt")
  
  if(file.exists(kmer_count_file) & (length(read_lines(kmer_count_file))!=0 )){
    
    # detect header line
    all_lines <- read_lines(kmer_count_file)
    hdr_idx <- which(startsWith(all_lines, "Background\tSignal\tSequence_length\tGap_position"))[1]
    cols <- strsplit(all_lines[hdr_idx], "\t")[[1]] #15
    
    # read the data rows after the header; keep exactly header’s columns
    df <- read_tsv(kmer_count_file, skip = hdr_idx, col_names = FALSE, show_col_types = FALSE) #18
    if(ncol(df)==18){
    names(df) <- c(cols[1:(length(cols)-1)], "local_max_type", "value1", cols[length(cols)], "value2")
    
    #remove rows from the end 
    
    df=df[-which(df$Background!=df$Background[1]),]
    
    df = df %>%
      separate(
        col = value2,
        into = c("value2", "local_hamming"),
        sep = "(?:\\t|\\\\t|\\s+)",  # tab, literal "\t", or whitespace
        extra = "merge",             # keep anything after the first split together
        fill  = "right"              # if missing flag, put NA
      ) %>%
      mutate(value2 = readr::parse_number(value2))
    
    }else{
      
      
      names(df) <- c(cols[1:(length(cols)-1)], "local_max_type", "value1", cols[length(cols)], "value2", "local_hamming")
      
      #remove rows from the end 
      
      df=df[-which(df$Background!=df$Background[1]),]
      
      
    }
    
    
    #test=df2 %>% filter(!is.na(local_hamming))
    
    #df2 has 16384 rows, seed length is 7, the number of possible k-mers is 4^7 16384
    
    df$Kmer %>% unique() %>% length() #16384
    
    #Load the seed realisations and their Hamming 1 distance k-mers
    kmers_Hamming1=read_delim(paste0(results_path, data$ID, ".tsv"), col_names=FALSE, delim="\t") #224
    names(kmers_Hamming1)=c("Kmer", "motif_match_score", "motif_match_rank")
    
    #Load the seed realisations and their Hamming 1&2 distance k-mers
    #kmers_Hamming12=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",
    #                                  "MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6", "_score_rank", 
    #                                  ".tsv"), col_names=FALSE) #1280
    #names(kmers_Hamming12)=c("Kmer", "score", "rank")
    
    #find . -maxdepth 1 -type f -printf '%s\t%f\n' | sort -n | numfmt --to=iec --field=1 --padding=7 | column -t
    
    
    match_ind=match(kmers_Hamming1$Kmer, df$Kmer)
    
    head(kmers_Hamming1$Kmer)
    head(df$Kmer[match_ind])
    
    kmer_counts=left_join(kmers_Hamming1, 
                          df %>% select("Kmer", "Signal count", "Background count", "Fold_change", "Local_max", "local_max_type", "local_hamming"), 
                          by = "Kmer")
    
    kmer_counts$`Corrected count`=kmer_counts$`Signal count`-(lambda*(N_sig/N_bg))*kmer_counts$`Background count`
    
    #Order based on the Corrected count and add rank
    
    kmer_counts <- kmer_counts %>%
      arrange(desc(`Corrected count`)) %>%
      mutate(`Corrected_count_rank` = row_number())
    
    write_delim(kmer_counts, #244
                file=paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_SELEX_counts.tsv"),
                delim="\t")
  
  } #SELEX data exists
  
  # Score based on jellyfish output 
  
  
  #kmer_signal_count_file=paste0("/scratch/project_2013895/SELEX/jellyfish/counts/", data$ID, "_signal_counts.txt")
  kmer_signal_count_file=paste0("/scratch/project_2013895/SELEX/jellyfish/filtered_counts/", data$ID, "_signal_counts.txt")
  
  
  if(file.exists(kmer_signal_count_file) ){
  
    kmer_background_count_file=paste0("/scratch/project_2013895/SELEX/jellyfish/filtered_counts/", data$ID, "_background_counts.txt")
    
    kmers_signal_counts=read_delim(kmer_signal_count_file, col_names=FALSE, delim=" ") #224
    kmers_background_counts=read_delim(kmer_background_count_file, col_names=FALSE, delim=" ") #224
    names(kmers_signal_counts)=c("Kmer","jellyfish_signal_count")
    names(kmers_background_counts)=c("Kmer","jellyfish_background_count")
    
    jellyfish_counts=left_join(kmers_signal_counts, kmers_background_counts, by="Kmer")
    rm(kmers_signal_counts, kmers_background_counts)
    
    kmers_Hamming1=read_delim(paste0(results_path, data$ID, ".tsv"), col_names=FALSE, delim="\t") #224
    names(kmers_Hamming1)=c("Kmer", "motif_match_score", "motif_match_rank")
    
    match_ind=match(kmers_Hamming1$Kmer, jellyfish_counts$Kmer)
    
    head(kmers_Hamming1$Kmer)
    head(jellyfish_counts$Kmer[match_ind])
    
    kmer_counts=left_join(kmers_Hamming1, 
                          jellyfish_counts, 
                          by = "Kmer")
    
    kmer_counts$`jellyfish Corrected count`=kmer_counts$jellyfish_signal_count-(lambda*(N_sig/N_bg))*kmer_counts$jellyfish_background_count
    
    #Order based on the Corrected count and add rank
    
    kmer_counts <- kmer_counts %>%
      arrange(desc(`jellyfish Corrected count`)) %>%
      mutate(`jellyfish Corrected_count_rank` = row_number())
    
    write_delim(kmer_counts, #244
                file=paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_jellyfish_counts.tsv"),
                delim="\t")
    
  } #jellyfish data exists
  
  #if both SELEX counts and jellyfish counts exist, create a combined table 
  
  if(file.exists(paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_jellyfish_counts.tsv")) &
     file.exists(paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_SELEX_counts.tsv"))){
    #kmer_counts_jellyfish=kmer_counts
    kmer_counts_jellyfish=read_delim(file=paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_jellyfish_counts.tsv"),
                                     delim="\t")
    
    kmer_counts_SELEX=read_delim(file=paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_SELEX_counts.tsv"),
                                 delim="\t")
    
    
    combined_counts=left_join(kmer_counts_SELEX, 
                              kmer_counts_jellyfish %>% select("Kmer",  "jellyfish_signal_count", "jellyfish_background_count", 
                                                               "jellyfish Corrected count","jellyfish Corrected_count_rank"), 
                              by = "Kmer")
    
    #combined_counts %>% select(`Signal count`, `Background count`, jellyfish_signal_count, jellyfish_background_count) %>% View()

    #isTRUE(all.equal(as.numeric(combined_counts$`Signal count`), as.numeric(combined_counts$jellyfish_signal_count)))
    #(combined_counts$`Signal count`==combined_counts$jellyfish_signal_count) %>% table()
    #identical((as.integer(combined_counts$`Background count`),as.integer(combined_counts$jellyfish_background_count))) %>% table()
    # tol <- 1e-8
    # all(abs(combined_counts$`Signal count` - combined_counts$jellyfish_signal_count) < tol, na.rm = TRUE)
    # x=combined_counts$`Signal count`
    # y=combined_counts$jellyfish_signal_count
    # bad <- which(!(abs(x - y) < tol) | xor(is.na(x), is.na(y)))
    # cbind(i = bad, x = x[bad], y = y[bad], diff = x[bad] - y[bad])
    #combined_counts[bad,] %>% View() These are palindromic k-mers? YES. spacek counts these twice
    
    #ind=head(which((combined_counts$`Signal count`==combined_counts$jellyfish_signal_count)))
    
    #combined_counts[ind,] %>% View()
    
    write_delim(combined_counts, #244
                file=paste0("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/",data$ID,"_combined_counts.tsv"),
                delim="\t")
  } #both count files exist
  }else{
    print("no lambda")
  } #lambda exists
}

  