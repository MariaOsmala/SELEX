library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")

args <- commandArgs(trailingOnly = TRUE)
arrays=as.numeric(args[1])
addition=as.numeric(args[2])
arrays=arrays+addition




metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

metadata=metadata %>% filter(!(experiment %in% "Methyl-HT-SELEX")) #3636

#failed_experiments=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/failed_experiments_Hamming1.tsv")
#failed_experiments_all=failed_experiments
#failed_experiments=failed_experiments %>% filter(State %in% c("CANCELLED", "FAILED"))

#failed_experiments=failed_experiments %>% filter(!(motif=="E2F8_Morgunova2015"))

#missing files

results=dir("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/")
results=gsub(".tsv", "", results) #3461

missing=metadata %>% filter(!(ID %in% results)) %>% pull(ID) #175
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

results_path="/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/" #3461
#3461+175=3636

#metadata=metadata %>% filter(ID %in% missing)

metadata=metadata %>% filter(!(ID=="E2F8_Morgunova2015"))



pseudocount=0.01

stream_canonical_and_hamming1_to_file <- function(iupac_seq, output_file) {
  iupac_map <- Biostrings::IUPAC_CODE_MAP
  bases <- c("A", "T", "C", "G")
  chars <- unlist(strsplit(iupac_seq, ""))
  base_options <- lapply(chars, function(ch) strsplit(iupac_map[[ch]], "")[[1]])
  
  seen <- new.env(hash = TRUE, parent = emptyenv())
  con <- file(output_file, open = "w")
  on.exit(close(con))
  
  # Function to get canonical form (smallest between k-mer and its reverse complement)
  canonical_form <- function(seq) {
    s <- DNAString(seq)
    min(as.character(s), as.character(reverseComplement(s)))
  }
  
  # Function to generate all Hamming distance 1 variants
  generate_hamming1 <- function(kmer) {
    chars <- unlist(strsplit(kmer, ""))
    variants <- character(0)
    
    for (i in seq_along(chars)) {
      for (alt in setdiff(bases, chars[i])) {
        mutated <- chars
        mutated[i] <- alt
        variant <- paste0(mutated, collapse = "")
        variants <- c(variants, canonical_form(variant))
      }
    }
    unique(variants)
  }
  
  # Recursive streaming and processing
  stream_recursive <- function(pos = 1, prefix = "") {
    if (pos > length(base_options)) {
      kmer <- prefix
      canon <- canonical_form(kmer)
      
      if (!exists(canon, envir = seen)) {
        assign(canon, TRUE, envir = seen)
        writeLines(canon, con)
      }
      
      # Now generate Hamming-1 variants
      for (v in generate_hamming1(kmer)) {
        if (!exists(v, envir = seen)) {
          assign(v, TRUE, envir = seen)
          writeLines(v, con)
        }
      }
      
    } else {
      for (base in base_options[[pos]]) {
        #base=base_options[[pos]][1]
        stream_recursive(pos + 1, paste0(prefix, base))
      }
    }
  }
  
  stream_recursive()
}


path_local_scratch=Sys.getenv("LOCAL_SCRATCH")
#dir(path_local_scratch)
#paste0("/scratch/project_2013895/SELEX/streamed_kmers/",metadata$ID[arrays],".tsv"))
stream_canonical_and_hamming1_to_file(metadata$seed[arrays], paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv"))



#iupac_seq="ATNGC"
#

#result=read_tsv(paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv"))
#result=read_tsv(paste0(path_local_scratch,"/",metadata$ID[arrays],"_score.tsv"), col_names = NULL)

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
  
  
index=arrays

#for(index in seq(start_ind, end_ind, 1)){ #0-9
  print(index)
  #index=which(metadata$ID=="MEIS1_DLX3_CAP-SELEX_TAAAGC40NGAA_AT_TGACANSNTAATTG_1_3")
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
  # paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv")
  kmer_file=paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv")
  output_file=paste0(path_local_scratch,"/",metadata$ID[arrays],"_score.tsv")
  
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
  
  score_kmers_with_motif(kmer_file=kmer_file, pssm=pssm, output_file=output_file) 
  
  #Do we get the same result: seems so
  
  # old=read.delim("/scratch/project_2013895/SELEX/kmers/MAX_HT-SELEX_TTCTGG40NTGC_KAE_RNCACGTGYN_1_3.tsv", header=FALSE)
  # new=read.delim("/scratch/project_2013895/SELEX/streamed_kmers/MAX_HT-SELEX_TTCTGG40NTGC_KAE_RNCACGTGYN_1_3_score.tsv", header=FALSE)
  # 
  # new=new[order(-new$V2), ] 
  # new$rank=1:nrow(new)
  # 
  # old==new  
  # 
  # (abs(old$V2- new$V2) < 1e-6) %>% table() #all TRUE
  
  #Sort the resulting file
  
  output_file_score_rank=paste0(path_local_scratch,"/",metadata$ID[arrays],".tsv")
  
  
  # Construct the shell command
  cmd <- sprintf("sort -k2,2nr %s | awk '{printf \"%%s\\t%%s\\t%%d\\n\", $1, $2, NR}' > %s", output_file, output_file_score_rank)
  
  # Run it
  system(cmd)
  system(paste0("rm ",output_file))
  
 
  
  system(paste0("cp ", output_file_score_rank," ",results_path ))
  
  
  
  #new=read.delim("/scratch/project_2013895/SELEX/streamed_kmers/MAX_HT-SELEX_TTCTGG40NTGC_KAE_RNCACGTGYN_1_3_score_rank.tsv", header=FALSE)
  