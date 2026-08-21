library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("motifmatchr")
library("TFBSTools")
library("rtracklayer")

args <- commandArgs(trailingOnly = TRUE)
arrays=as.numeric(args[1])

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

start_ind=arrays*10+1 #100
end_ind=(arrays+1)*10 #100
len=10 #100

#3933 motifs
if(end_ind>393){ # 398 
  end_ind=3933
}

# This code could also score the k-mers based on PWM
# sort the k-mers based on the scores and give rank 

# Function to expand IUPAC ----------------------------------------------

iupac_map <- Biostrings::IUPAC_CODE_MAP
expand_iupac <- function(seq) {
  chars <- strsplit(as.character(seq), "")[[1]]
  bases <- lapply(chars, function(ch) strsplit(iupac_map[[ch]], "")[[1]])
  combos <- expand.grid(bases, stringsAsFactors = FALSE)
  apply(combos, 1, paste0, collapse = "")
}

# Function to generate Hamming distance 1 variants for each k-mer ----------------------------------------------

generate_hd1 <- function(seq) {
  bases <- c("A", "C", "G", "T")
  chars <- strsplit(seq, "")[[1]]
  variants <- c()
  
  for (i in seq_along(chars)) {
    #i=1
    for (b in bases) {
      #b=bases[2]
      if (b != chars[i]) { #Do not convert the base to itself
        temp <- chars
        temp[i] <- b
        variants <- c(variants, paste(temp, collapse = ""))
      }
    }
  }
  return(unique(variants))
}


# Function to create PWM from PFM -----------------------------------------


PPMp <- function(C, p){
  (C + p / nrow(C)) / matrix( (colSums(C) + p),nrow=4, ncol=ncol(C),byrow=TRUE)
  
}


# Position-specific scoring matrix ----------------------------------------

#Position specific scoring matrix
PSSM <- function(C, B) log2(C / B) #log

score_sequence <- function(seq, pssm){
  #seq="CAGTGTGGTCGC"
  seq=sapply(as.character(seq), strsplit, "")[[1]]
  sum(as.numeric(mapply(function(base, pos) pssm[base, pos], seq, seq_along(seq))))
}




for(index in seq(start_ind, end_ind, 1)){ #0-9
  
  data=metadata[index,] 
  
  seed=data$seed
  
  
  # Create an IUPAC DNA string ----------------------------------------------
  seq <- DNAString(seed)
  
  
  
  resolved_seqs <- expand_iupac(seq)
  
  resolved_seqs <- unique(resolved_seqs)
  
  
  
  all_kmers <- unique(unlist(lapply(resolved_seqs, generate_hd1)))
  
  
  
  ## Include original k-mers to the list -------------------------------------
  all_kmers <- unique(c(resolved_seqs, all_kmers)) #4864
  
  
  # Check that the hamming distances are all 1 ------------------------------
  
  variants_list=lapply(resolved_seqs, generate_hd1)
  
  hm=list()
  for(i in 1:length(resolved_seqs)){
    hm[[i]]=stringdist(resolved_seqs[i], variants_list[[i]], method = "hamming")
  }
  
  unique(sapply(hm, unique)) #1 YES THEY ARE!
  
  
  
  # Remove reverse complements ----------------------------------------------
  
  
  # Step 1: Convert to DNAStringSet
  kmer_set <- DNAStringSet(all_kmers)
  
  # Step 2: Compute reverse complements
  rc_set <- reverseComplement(kmer_set)
  
  # Step 3: Convert both to character vectors
  kmer_chars <- as.character(kmer_set)
  rc_chars <- as.character(rc_set)
  
  # Step 4: Define canonical form — lexicographically smaller of the pair
  canonical <- pmin(kmer_chars, rc_chars)
  
  table(canonical) %>% table()
  #1(without reverse complements)    2 (with reverse complements) 
  #16 2424 
  
  # Step 5: Keep only unique canonical forms
  unique_kmers <- unique(canonical)
  
  # Result:
  length(unique_kmers)  # Number of unique k-mers after collapsing rev-comp pairs 2440
  #head(unique_kmers)
  
  
  # Score the k-mers by pwm -------------------------------------------------
  
  pcm=as.matrix(read.table(file=gsub("../../","/projappl/project_2006203/TFBS/",data$filename)))
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  pfm <- TFBSTools::PFMatrix(ID=data$ID, 
                             strand="+", 
                             bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                             profileMatrix=pcm
  )
  
  pseudocount=0.01
  
  
  #Check that two approaches gives the same result
  pwm=TFBSTools::toPWM(pfm, type="prob", pseudocounts=pseudocount, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  PWM=PPMp(pfm@profileMatrix,pseudocount)
  
  #NOTICE THAT THESE ARE LOG2 PSSM matrices
  pssm_test=TFBSTools::toPWM(pfm, type="log2probratio", 
                             pseudocounts=pseudocount, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  pssm=PSSM(PWM, 0.25)
  #Now these are the same
  #pssm_test@profileMatrix==pssm
  
  # I need to score both the k-mer and its reverse complement and take the best score
  
  # Step 1: Convert to DNAStringSet
  kmer_set <- DNAStringSet(unique_kmers)
  
  # Step 2: Compute reverse complements
  rc_set <- reverseComplement(kmer_set)
  
  # Step 3: Convert both to character vectors
  kmer_chars <- as.character(kmer_set)
  rc_chars <- as.character(rc_set)
  
  scores=data.frame(f_scores=sapply(kmer_chars,function(x) score_sequence(x, pssm_test@profileMatrix)),
                    rc_scores=sapply(rc_chars,function(x) score_sequence(x, pssm_test@profileMatrix)))
  
  scores_max=apply(scores, 1, max)
  
  result=data.frame(V1=unique_kmers, score=scores_max)
  
  #order 
  result=result[order(-result$score), ] 
  #add rank
  result$rank=1:nrow(result)  
  
  #Another function to score
  #score <- function(S, seq) {
  #  
  #  sum(S[cbind(match(strsplit(as.character(seq), "")[[1]], rownames(S)), 1:ncol(S))])
  #}
  
  #score(pssm_test@profileMatrix, reverseComplement(DNAStringSet("CAATTAACCTGTCA"))) # 21.8616
  
  #What is the pseudocount used here
  universalmotif::scan_sequences(
    universalmotif::convert_motifs(PWM), 
    reverseComplement(DNAStringSet("CAATTAACCTGTCA")))$score #21.853
  #moods gives 10.71408, this is likely the log2 score
  universalmotif::scan_sequences(
    universalmotif::convert_motifs(pwm), 
    reverseComplement(DNAStringSet("CAATTAACCTGTCA")))$score #10.709 THIS IS WRONG
  universalmotif::scan_sequences(
    universalmotif::convert_motifs(pfm), 
    reverseComplement(DNAStringSet("CAATTAACCTGTCA")))$score #21.852
  
  universalmotif::scan_sequences(
    universalmotif::convert_motifs(pssm_test), 
    reverseComplement(DNAStringSet("CAATTAACCTGTCA")))$score #21.853
  
  test=Biostrings::matchPWM(PWM, reverseComplement(DNAString("CAATTAACCTGTCA")), with.score=TRUE )
  test@elementMetadata$score #10.71408 This is wrong
  
  test=Biostrings::matchPWM(pssm, reverseComplement(DNAString("CAATTAACCTGTCA")), with.score=TRUE )
  test@elementMetadata$score #21.8616 This is correct
  
  
  
  #What are the matche scores in the genome 
  #MOODS=rtracklayer::import("/scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct_processed/MOODS_bigbed/MEIS1_DLX3_CAP-SELEX_TAAAGC40NGAA_AT_TGACANSNTAATTG_1_3_top.bed")
  #Maximum moods score 15.1533
  #max(MOODS$score)
  #Convert to log2 score
  #max(MOODS$score)/log(2) #21.86161
  
  
  # Do we get the same score by MOODS ---------------------------------------
  
  #There can be 0 1 or 2 matches
  
  #sapply(matchMotifs(pwm, DNAStringSet(kmer_chars), out="positions")[[1]], length) %>% table()#All 1
  
  #matches=as.data.frame(unlist(matchMotifs(pwm, DNAStringSet(result$V1), #pwm gives wrong scores
  #                                         p.cutoff = 5e-04,
  #                                         out="positions")[[1]]))
  
  #Do not know what is the pseudocount used here, maybe 1
  matches=as.data.frame(unlist(matchMotifs(pfm, DNAStringSet(result$V1), #pwm gives wrong scores
                                           p.cutoff = 5e-03,bg="even",
                                           out="positions")[[1]]))
  
  nrow(matches)==nrow(result)
  #matches are log2 scores
  #matches$score=matches$score/log2(exp(1))
  matches$seq=result$V1
  
  matches=matches[order(-matches$score),]
  
  
  
  write_delim(data.frame(V1=unique_kmers), 
              file=paste0("/scratch/project_2013895/SELEX/kmers/",
                          data$ID, ".tsv"),
              col_names = FALSE, delim="\t")
}
