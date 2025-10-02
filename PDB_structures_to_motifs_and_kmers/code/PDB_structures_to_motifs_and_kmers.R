library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

rm(list=ls())


PDB <- read_excel("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/PDB_structures_160725_PDB_link.xlsx", 
                                               skip = 10)
homodimers = PDB %>% filter(Type=="homodimer")
heterodimers = PDB %>% filter(Type=="heterodimer")
heterotrimers = PDB %>% filter(Type=="heterotrimer")

homodimers=homodimers%>%
  mutate(`DNA sequence` = strsplit(`DNA sequence`, ";\\s*")) %>%  # split into list column
  unnest(`DNA sequence`)

metadata <- read_tsv("/projappl/project_2013895/motif_metadata/metadata_final.tsv", col_names=TRUE)

metadata_monomers=metadata %>% filter(experiment %in% c("HT-SELEX", "Methyl-HT-SELEX"))

(metadata_monomers$symbol %in% homodimers$TFs) %>% table()
#FALSE  TRUE 
#1924   111 

(homodimers$TFs %in% metadata_monomers$symbol) %>% table()
#FALSE  TRUE 
#59    44 

homodimers=homodimers %>% filter(TFs %in% metadata_monomers$symbol)

TFs=names(homodimers$TFs %>% table() %>%sort(decreasing=TRUE))

homodimers=homodimers %>% mutate(motif="") %>% relocate(motif, .after=`DBD type`)
homodimers=homodimers %>% mutate(motif_type="") %>% relocate(motif_type, .after=`motif`)
homodimers=homodimers %>% mutate(seed="") %>% relocate(seed, .after=`motif_type`)
homodimers=homodimers %>% mutate(motif_length=NA) %>% relocate(motif_length, .after=`seed`)
homodimers=homodimers %>% mutate(kmer="") %>% relocate(kmer, .after=`motif_length`)
homodimers=homodimers %>% mutate(orientation="") %>% relocate(orientation, .after=`kmer`)
homodimers=homodimers %>% mutate(scaled_score=NA) %>% relocate(scaled_score, .after=`orientation`)
homodimers=homodimers %>% mutate(rank=NA) %>% relocate(rank, .after=`scaled_score`)
homodimers=homodimers %>% mutate(orig_score=NA) %>% relocate(orig_score, .after=`rank`)
homodimers=homodimers %>% mutate(kmer_alignment_start=NA) %>% relocate(kmer_alignment_start, .after=`orig_score`)
homodimers=homodimers %>% mutate(kmer_alignment_end=NA) %>% relocate(kmer_alignment_end, .after=`kmer_alignment_start`)
homodimers=homodimers %>% mutate(structure_sequence_alignment_start=NA) %>% relocate(structure_sequence_alignment_start, .after=`kmer_alignment_end`)
homodimers=homodimers %>% mutate(structure_sequence_alignment_end=NA) %>% relocate(structure_sequence_alignment_end, .after=`structure_sequence_alignment_start`)
homodimers=homodimers %>% mutate(alignment_score=NA) %>% relocate(alignment_score, .after=`structure_sequence_alignment_end`)

library(Biostrings)

for(TF in TFs){
  #TF=TFs[1]
  
  
  structure_sequences=homodimers %>% filter(TFs==TF) %>% pull(`DNA sequence`) %>% unique() %>% DNAStringSet()
  
  motifs=metadata_monomers %>% filter(symbol==TF)
  print(paste0(TF, ": ",nrow(motifs), " motifs"))
  motif_consensus=motifs$consensus
  motif_seed=motifs$seed
  
  source("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/code/functions.R")
  
  #Need to go through all motifs
  
  results_list=list()
  
  for(motif in motifs$ID){
    print(motif)
    #motif=motifs$ID[1]
    # Extract k-mers 
    
    #3263
    #kmers_scores=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",motifs$ID,".tsv"), col_names = FALSE)
    
    #50912
   kmers_scores=try(read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",motif,"_score_rank.tsv"), col_names = FALSE),
        silent=TRUE)
   if(!("try-error" %in% class(kmers_scores))){
      kmers_scores_orig=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/",motif,"_score_rank.tsv"), col_names = FALSE)
   }else{
     kmers_scores=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",motif,".tsv"), col_names = FALSE)
     kmers_scores_orig=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/",motif,".tsv"), col_names = FALSE)
   }
    
    names(kmers_scores)=c("kmers", "scores", "rank")
    names(kmers_scores_orig)=c("kmers", "scores", "rank")
    kmers_scores=kmers_scores %>% mutate(orig_scores=kmers_scores_orig$scores) %>% relocate(orig_scores, .after = scores)
    rm(kmers_scores_orig)
    
    nchar(kmers_scores$kmers) %>% unique() #10
    
    kmers=DNAStringSet(kmers_scores$kmers)
    
    pdict=c(kmers, reverseComplement(kmers))
    
    structure_sequence_list=list()
    
    for(i in 1:length(structure_sequences)){
      #i=1
      subject=structure_sequences[[i]]
      print(subject)
      match_ind=whichPDict(pdict, subject,
                 max.mismatch=0, min.mismatch=0, with.indels=TRUE, fixed=TRUE,
                 algorithm="auto", verbose=FALSE)
      if(length(match_ind)!=0){
        
        if(length(match_ind)>1){
          print("several match_inds!")
        }
        
        if( match_ind <= nrow(kmers_scores) ){
           match_with_orig=TRUE
           alignment=align_with_revcomp(pdict[match_ind],subject)$best
           kmer=pdict[match_ind]
           
           
        }else{
          match_with_orig=FALSE #Match with reverse complement
          match_ind=match_ind-nrow(kmers_scores)
          alignment=align_with_revcomp(pdict[match_ind],subject)$best
          kmer=pdict[match_ind]
        }
        
        structure_sequence_list[[as.character(subject)]]=list(match_ind=match_ind,
          match_with_orig=match_with_orig,
                                   alignment=alignment,
                                    alignment_score=alignment@score,
                                   kmer=kmer
          
        )
        
        
        
      }else{
        
        structure_sequence_list[[as.character(subject)]]=list(match_ind=NA,
          match_with_orig=NA,
                                   alignment=NA,
          alignment_score=NA,
                                   kmer=NA)
        
      }
      
      
    }
    
    results_list[[motif]]=structure_sequence_list
  
  } #over motifs

  #Which sequences have matches to the motif k-mers
  
  results=as.data.frame(do.call(rbind, lapply(results_list, function(x) sapply(x, function(y) !is.na(y$match_ind)))))
  
  motif_matches=results %>% filter(CTAAACGGGCAATTAG==TRUE) %>% row.names()
  #What are the scores, ranks and alignment scores
  
  lapply( results_list[[motif_matches]], function(x) x[["CTAAACGGGCAATTAG"]])
  
  results %>% filter(CTAATTGCTACCGTTTAG==TRUE) %>% row.names()
  
  
  
  ind=which(homodimers$TFs==TF & homodimers$`DNA sequence`==as.character(subject))
  
  homodimers$motif[ind]=motifs$ID
  homodimers$motif_type[ind]=motifs$type
  homodimers$seed[ind]=motifs$seed
  homodimers$motif_length[ind]=motifs$length
  homodimers$kmer[ind]=as.character(kmer)
  homodimers$kmer_orig[ind]=match_with_orig
  homodimers$scaled_score[ind]=kmers_scores$scores[match_ind]
  homodimers$rank[ind]=kmers_scores$rank[match_ind]
  homodimers$orig_score[ind]=kmers_scores$orig_scores[match_ind]
  homodimers$kmer_alignment_start[ind]=start(alignment@pattern@range) #1
  homodimers$kmer_alignment_end[ind]=end(alignment@pattern@range) #15
  homodimers$structure_sequence_alignment_start[ind]= start(alignment@subject@range) #4
  homodimers$structure_sequence_alignment_end[ind]=end(alignment@subject@range)  #18
  homodimers$alignment_score[ind]=alignment@score

  
  
  library(Biostrings)
  library(dplyr)
  library(readr)
  # assumes you already have: TFs, homodimers, metadata_monomers, align_with_revcomp()
  
  all_hits <- list()
  
  #RFX1 as uracil in the sequence
  
  for (TF in TFs[29:length(TFs)]) {
    #TF=TFs[2]
    structure_sequences <- homodimers %>%
      filter(TFs == TF) %>%
      pull(`DNA sequence`) %>%
      unique() %>%
      DNAStringSet()
    
    motifs <- metadata_monomers %>% filter(symbol == TF)
    motifs<- motifs %>% filter(experiment=="HT-SELEX")
    message(paste0(TF, ": ", nrow(motifs), " motifs"))
    
    hits <- list()  # collect rows for this TF
    
    for (motif in motifs$ID) {
      #motif=motifs$ID[1]
      message(motif)
      
      # Try fast (50912) source; fall back to dense set
      kmers_scores <- try(
        read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/", motif, "_score_rank.tsv"),
                   col_names = FALSE),
        silent = TRUE
      )
      if (!inherits(kmers_scores, "try-error")) {
        kmers_scores_orig <- read_delim(
          paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/", motif, "_score_rank.tsv"),
          col_names = FALSE
        )
      } else {
        kmers_scores <- read_delim(
          paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/", motif, ".tsv"),
          col_names = FALSE
        )
        kmers_scores_orig <- read_delim(
          paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/", motif, ".tsv"),
          col_names = FALSE
        )
      }
      
      names(kmers_scores) <- c("kmers", "score_scaled", "rank")
      names(kmers_scores_orig) <- c("kmers", "score_orig", "rank")
      kmers_scores <- kmers_scores %>%
        mutate(score_orig = kmers_scores_orig$score_orig) %>%
        relocate(score_orig, .after = score_scaled)
      rm(kmers_scores_orig)
      
      #nchar(kmers_scores$kmers) %>% table()
      
      kmers <- DNAStringSet(kmers_scores$kmers)
      pdict <- c(kmers, reverseComplement(kmers))   # forward + RC dictionary
      
      for (i in seq_along(structure_sequences)) {
        #i=1
        subject <- structure_sequences[[i]]
        
        # IUPAC-aware exact match search (no indels here)
        match_ind <- whichPDict(
          pdict, subject,
          max.mismatch = 0, min.mismatch = 0,
          with.indels = FALSE, fixed = FALSE, algorithm = "auto"
        )
        
        if (length(match_ind) == 0L) next
        
        if (length(match_ind) > 1L) {
          # if multiple patterns hit, take the first (or handle all—up to you)
          match_ind <- match_ind[1]
        }
        
        # Determine orientation and map to original forward k-mer index
        if (match_ind <= nrow(kmers_scores)) {
          match_with_orig <- TRUE
          idx <- match_ind
          orientation <- "forward"
        } else {
          match_with_orig <- FALSE
          idx <- match_ind - nrow(kmers_scores)
          orientation <- "reverse_complement"
        }
        
        # Best gapped alignment (IUPAC-aware if your align_with_revcomp uses the IUPAC matrix)
        aln <- align_with_revcomp(kmers[idx], subject)$best
        
        # Collect a tidy row
        hits[[length(hits) + 1L]] <- tibble(
          TF = TF,
          motif_id = motif,
          motif_length=metadata_monomers %>% filter(ID==motif) %>% pull(length),
          seed=metadata_monomers %>% filter(ID==motif) %>% pull(seed),
          motif_type=metadata_monomers %>% filter(ID==motif) %>% pull(type),
          structure_sequence = as.character(subject),
          kmer = as.character(kmers[idx]),
          orientation = orientation,
          kmer_rank = kmers_scores$rank[idx],
          kmer_score_scaled = kmers_scores$score_scaled[idx],
          kmer_score_orig = kmers_scores$score_orig[idx],
          alignment_score = score(aln),
          kmer_alignment_start=start(alignment@pattern@range), #1
          kmer_alignment_end=end(alignment@pattern@range), #15
          structure_sequence_alignment_start= start(alignment@subject@range), #4
          structure_sequence_alignment_end=end(alignment@subject@range)  #18
        )
      }
    } #motifs
    
    all_hits[[TF]] <- if (length(hits)) bind_rows(hits) else tibble()
  
 
  if(nrow(all_hits[[TF]])>0){
    df=all_hits[[TF]]  %>% group_by(structure_sequence) %>% nest()
    
    for(j in 1:length(df$structure_sequence)){
      #j=1
      seq=df$structure_sequence[j]
      max_ind=which(df$data[[j]]$alignment_score==max(df$data[[j]]$alignment_score))
      
      df$data[[j]]=df$data[[j]][max_ind,]
      min_ind=which(df$data[[j]]$kmer_rank==min(df$data[[j]]$kmer_rank))
      df$data[[j]]=df$data[[j]][min_ind,]
      
      if(nrow(df$data[[j]])>1){
        print("multiple matches")
      }
      
      ind=which(homodimers$TFs==TF & homodimers$`DNA sequence`==seq)
      
      homodimers$motif[ind]=df$data[[j]]$motif_id
      homodimers$motif_type[ind]=df$data[[j]]$motif_type
      homodimers$seed[ind]=df$data[[j]]$seed
      homodimers$motif_length[ind]=df$data[[j]]$motif_length
      homodimers$kmer[ind]=df$data[[j]]$kmer
      homodimers$orientation[ind]=df$data[[j]]$orientation
      homodimers$scaled_score[ind]=df$data[[j]]$kmer_score_scaled
      homodimers$rank[ind]=df$data[[j]]$kmer_rank
      homodimers$orig_score[ind]=df$data[[j]]$kmer_score_orig
      homodimers$kmer_alignment_start[ind]=df$data[[j]]$kmer_alignment_start
      homodimers$kmer_alignment_end[ind]=df$data[[j]]$kmer_alignment_end
      homodimers$structure_sequence_alignment_start[ind]= df$data[[j]]$structure_sequence_alignment_start
      homodimers$structure_sequence_alignment_end[ind]=df$data[[j]]$structure_sequence_alignment_end
      homodimers$alignment_score[ind]=df$data[[j]]$alignment_score
    }
  
  }
  
}

write_delim(homodimers,file="/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/homodimers.tsv", delim="\t", col_names=TRUE)  
    
  
  
  
  