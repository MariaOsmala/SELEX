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
arrays=as.numeric(args[1]) #1-360
addition=as.numeric(args[2])
arrays=arrays+addition

#length=10 
length=10 #starting from 231 2301-2310
start_ind=(arrays-1)*length+1 #100
end_ind=arrays*length #100

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635

#which(metadata$seed=="NNNACGANNNNNNTCGTNNN") #992
#which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3")



#sort based on the seed length

metadata=metadata[order(metadata$length),]


#results_path="/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/" #
results_path="/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds_June2026/" #3595, one missing?

#Filter out those which do not have either SELEX signal or background
metadata=metadata %>% filter(!(is.na(CSC_SELEX_filename) | is.na(CSC_SELEX_background_filename))) #3596

#which(metadata$ID=="TGIF2_TBX21_TGGGTC40NAAG_YYI_NAGGTGTCAWN_1_3") 682
#which(metadata$ID=="Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3") #373
#which(metadata$ID=="BARHL2_HT-SELEX_TCCAGT40NGAC_AI_NNTAATTGNN_1_3") #377
#which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3") #2333
#which(metadata$ID=="MEIS2_HT-SELEX_TCGCTG20NGA_AA_TTGACAGSTGTCAA_2_4") #1442
#which(metadata$ID=="MEIS1_DLX3_CAP-SELEX_TAAAGC40NGAA_AT_TGACANSNTAATTG_1_3")
#which(metadata$ID=="MEIS1_HOXB13_CAP-SELEX_TAAGGA40NCCG_AAB_SYMRTAAANCTGTCA_1_3")
#which(metadata$ID=="FOXO1_ETV4_CAP-SELEX_TCCGTA40NGCC_AS_RWMAACAGGAARNN_1_3b0")
#which(metadata$ID=="JDP2_HT-SELEX_TGTTCA20NGA_AE_ATGASTCAT_1_3") #216

# which(metadata$ID=="GCM2_TBX21_CAP-SELEX_TGTTCA40NCTG_AAB_NNNCGGGNNNGGTGTNN_1_3") #2601
# which(metadata$ID=="UBP1_HT-SELEX_TGGCCT40NTAT_KT_NNCCGGNNNNNNCCGGNN_2_4")        #2993    
# which(metadata$ID=="ETS2_HT-SELEX_TGTGAT40NGCA_KX_NACCGGANNNNNNTCCGGTN_1_4")      #3204
# which(metadata$ID=="ETV3_HT-SELEX_TTATTC40NGAA_KV_NNAGGAANNNNNNNTTCCTNN_1_3")     #3393
# which(metadata$ID=="HOXB2_ETV7_CAP-SELEX_TGAGGT40NTCC_AY_TAATKNNNNGNNNNNNCTTCCNN_1_3") #3473    
# which(metadata$ID=="JDP2_EVX1_TGCTTT40NAGA_YPIIII_TGACTCNNNNNNNNNNNTAATTA_1_3")   #3495    
# which(metadata$ID=="POU4F3_CREM_TAAGTT40NGAT_YLIII_GACGTCNNNNNNNNNNNNTATGCA_1_3") #3522   

if(end_ind > nrow(metadata)){
  end_ind=nrow(metadata)
}

#metadata$index=1:nrow(metadata)
#results=dir("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/")
#results_signal=results[grep("signal", results)] #3474
#results_background=results[grep("background", results)] #3573
#missing_signal=metadata %>% filter(!(ID %in% gsub("_signal.txt", "", results_signal))) #820
#missing_background=metadata %>% filter(!(ID %in% gsub("_background.txt", "", results_background))) #897
#missing_background=metadata %>% filter(!(ID %in% results_background))
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
# ls | wc -l #3596

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

pseudocount=0.01

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
  #seq=as.character(kmer_dna[[1]])
  seq=sapply(as.character(seq), strsplit, "")[[1]]
  sum(as.numeric(mapply(function(base, pos) pssm[base, pos], seq, seq_along(seq))))
}

score_pssm_expectedN <- function(pssm, seq, bg = c(A=.25, C=.25, G=.25, T=.25)) {
  #seq=as.character(kmer_dna[[1]])
  seq=sapply(as.character(seq), strsplit, "")[[1]]
  stopifnot(length(seq) == ncol(pssm))
  stopifnot(all(names(bg) %in% rownames(pssm)))
  bg <- bg[rownames(pssm)] / sum(bg)
  
  # per-position expected score under bg
  exp_col <- as.numeric(t(bg) %*% pssm[rownames(pssm), , drop = FALSE])
  
  s <- 0
  for (i in seq_along(seq)) {
    b <- seq[i]
    if (b %in% rownames(pssm)) {
      s <- s + pssm[b, i]
    } else if (b == "N") {
      s <- s + exp_col[i]
    } else {
      stop("Unexpected base: ", b)
    }
  }
  as.numeric(s)
}


#tmp=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv") #3573
#which(tmp$ID==data$ID)

#index=2333

# which(metadata$ID=="GCM2_TBX21_CAP-SELEX_TGTTCA40NCTG_AAB_NNNCGGGNNNGGTGTNN_1_3") #2601
# which(metadata$ID=="UBP1_HT-SELEX_TGGCCT40NTAT_KT_NNCCGGNNNNNNCCGGNN_2_4")        #2993    
# which(metadata$ID=="ETS2_HT-SELEX_TGTGAT40NGCA_KX_NACCGGANNNNNNTCCGGTN_1_4")      #3204
# which(metadata$ID=="ETV3_HT-SELEX_TTATTC40NGAA_KV_NNAGGAANNNNNNNTTCCTNN_1_3")     #3393
# which(metadata$ID=="HOXB2_ETV7_CAP-SELEX_TGAGGT40NTCC_AY_TAATKNNNNGNNNNNNCTTCCNN_1_3") #3473    
# which(metadata$ID=="JDP2_EVX1_TGCTTT40NAGA_YPIIII_TGACTCNNNNNNNNNNNTAATTA_1_3")   #3495    
# which(metadata$ID=="POU4F3_CREM_TAAGTT40NGAT_YLIII_GACGTCNNNNNNNNNNNNTATGCA_1_3") #3522  

for(index in seq(start_ind, end_ind, 1)){ #0-9

  
  print(index)
  data=metadata[index,] 
  seed=data$seed
  
  # TNS CCN NNG GSN NNNNNNNNN NNN NCA CGTGN 6*3+5+9 32
  # NCA CGT GNN NNN NNNNNNNNN CCC NNN GGCNA
  # NCA CGT GNN NNN NNNNNNNNN CCC NNN GGCNA
  res <- try(parse_spacek_report(paste0("/scratch/project_2013895/SELEX/spacek/lambda/", data$ID, ".txt")), silent=TRUE)
  
  if(class(res)!="try-error"){
    if(!is.na(res$lambda)){
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
        mutate(`Corrected_count_rank` = row_number()) %>% select(Kmer, `Corrected_count`, `Corrected_count_scaled`, Corrected_count_rank)
      
      #score k-mers with motif
      
      pcm=as.matrix(read.table(file=gsub("../../","/projappl/project_2006203/TFBS/",data$filename))) #31 columns
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
    
    
      kmer_dna <- sapply(as.list(counts$Kmer), DNAString)
      rc_kmer_dna <- sapply(kmer_dna, function(x) reverseComplement(DNAString(x)))
          
      # I need to score both the k-mer and its reverse complement and take the best score
      # Write kmer and score
      
      #TFAP2C_MAX_CAP-SELEX_TTAGTC40NTCC_AY_TNSCCNNNGGSNNNNNNNNNNNNNNCACGTGN_1_3"
      #TFAP2C_MAX_CAP-SELEX_TTAGTC40NTCC_AY_TNSCCNNNGGSNNNNNNNNNNNNNNCACGTGN_1_3.pfm"
      #                                     TNSCCNNNGG SNNNNNNNNN NNNNNCACGT GN" #32
      
     kmer_score=data.frame(orig=rep(NA,length(kmer_dna)), rc=rep(NA,length(kmer_dna)))
     kmer_score$orig=sapply(kmer_dna, function(x) score_pssm_expectedN(pssm, as.character(x)))
     kmer_score$rc=sapply(rc_kmer_dna, function(x) score_pssm_expectedN(pssm, as.character(x)))
    
     counts$motif_match_score=apply(kmer_score, 1, max)
     
     #max(counts$motif_match_score) can be negative
     counts$motif_match_score_scaled=counts$motif_match_score/abs(max(counts$motif_match_score))
     
     
     counts <- counts %>%
       arrange(desc(`motif_match_score_scaled`)) %>%
       mutate(`motif_match_score_rank` = row_number()) %>% select(Kmer, Corrected_count, Corrected_count_scaled, Corrected_count_rank, motif_match_score, motif_match_score_scaled, motif_match_score_rank)
     
     counts <- counts %>%
       arrange(`Corrected_count_rank`) 
     
     write_delim(counts, #244
                 file=paste0(results_path,data$ID,".tsv"),
                 delim="\t")
   

    }
  }else{
    print("no lambda")
  } #lambda exists
}

  