library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(Biostrings)
library(dplyr)
library(readr)

rm(list=ls())


#PDB <- read_excel("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/PDB_structures_160725_PDB_link.xlsx", 
PDB <- read_excel("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/PDB_structures_with_info.xlsx") 
                                               
homodimers = PDB %>% filter(Type=="homodimer")
heterodimers = PDB %>% filter(Type=="heterodimer")
heterotrimers = PDB %>% filter(Type=="heterotrimer")


test=homodimers %>% select(TFs, TFs_alternative,   uniprot1_hgnc_symbol_external_gene_name,uniprot2_hgnc_symbol_external_gene_name,  
                           uniprot3_hgnc_symbol_external_gene_name,      uniprot4_hgnc_symbol_external_gene_name,
                            `DNA sequence`, DNA1_sequence, DNA2_sequence,DNA3_sequence,DNA4_sequence,
                     
                      nakb_DNA1_description, nakb_DNA2_description, nakb_DNA3_description, nakb_DNA4_description,
                      #DNA1_molecule, DNA2_molecule,  DNA3_molecule,  DNA4_molecule,           
                      DNA1_organism, DNA2_organism, DNA3_organism, DNA4_organism,
                      nakb_DNA1_source,          nakb_DNA2_source,          nakb_DNA3_source,          nakb_DNA4_source,
                      DNA1_chains, DNA2_chains,  DNA3_chains,               DNA4_chains,
                      nakb_DNA1_chains,       nakb_DNA2_chains,          nakb_DNA3_chains,          nakb_DNA4_chains     
                      #struct_title, struct_pdbx_keywords, struct_text_keywords,  NAKBTitle,  NAKBna,  NAKBprot, 
                     
                      )                 

is.na(test) %>% table()

homodimers <- homodimers %>%
  mutate(
    across(where(is.character), ~ replace_na(.x, "")),
    across(where(is.factor),    ~ forcats::fct_explicit_na(.x, na_level = ""))
  )



#The sequences in the following columns seems to match
seq_cols  <- paste0("DNA", 1:4, "_sequence")
desc_cols <- paste0("nakb_DNA", 1:4, "_description")

# helper: normalize strings for matching
norm <- function(x) {
  x |>
    str_remove("\\(\\s*\\d+\\s*-?\\s*MER\\s*\\)") |>  # drop "(11-MER)" etc
    str_replace_all("\\s+", "") |>                    # remove spaces/tabs
    toupper()
}

df_check <- homodimers %>%
  # normalize sequences and descriptions (keep originals untouched)
  mutate(
    across(all_of(seq_cols),  norm, .names = "{.col}__c"),
    across(all_of(desc_cols), norm, .names = "{.col}__c")
  ) %>%
  rowwise() %>%
  mutate(
    # make one blob with all descriptions in this row
    .desc_blob = paste(na.omit(c_across(all_of(paste0(desc_cols, "__c")))), collapse = " "),
    # for each DNA*_sequence: does it occur in any description?
    across(all_of(paste0(seq_cols, "__c")),
           ~ ifelse(is.na(.x) || .x == "", NA, str_detect(.desc_blob, fixed(.x))),
           .names = "{.col}_in_desc")
  ) %>%
  ungroup() %>%
  # clean helper columns and names
  rename_with(~ sub("__c", "", .x), ends_with("_in_desc")) %>%
  select(-ends_with("__c"), -.desc_blob) %>%
  # row-level summary: all present? and which sequences are missing?
  rowwise() %>%
  mutate(
    all_dna_sequences_matched =
      if_all(ends_with("_in_desc"), ~ replace_na(.x, TRUE)),
    missing_dna_sequences = {
      seq_vals <- c_across(all_of(seq_cols))
      hit_vals <- c_across(ends_with("_in_desc"))
      miss <- seq_vals[!is.na(hit_vals) & hit_vals == FALSE]
      if (length(miss)) paste(miss, collapse = "; ") else NA_character_
    }
  ) %>%
  ungroup() %>% select( all_of(c(seq_cols, desc_cols,c("DNA1_sequence_in_desc", "DNA2_sequence_in_desc",                  
                         "DNA3_sequence_in_desc", "DNA4_sequence_in_desc", "all_dna_sequences_matched",              
                        "missing_dna_sequences"))))


#are the sequences reverse complments of each other 
library(Biostrings)



is_revcomp <- function(a, b) {
#  print(a)
#  print(b)
  A <- DNAStringSet(a)
  B <- DNAStringSet(b)
  as.character(reverseComplement(A)) == as.character(B)
}


DNAStringSet(homodimers$DNA1_sequence)

is_atcg1 <- !(homodimers$DNA1_sequence=="") & str_detect(homodimers$DNA1_sequence, regex("^[ACGT]+$", ignore_case = TRUE))
is_atcg2 <- !(homodimers$DNA2_sequence=="") & str_detect(homodimers$DNA2_sequence, regex("^[ACGT]+$", ignore_case = TRUE))

which(is_atcg1 & is_atcg2)

# Example with columns
DNA1_is_rc_of_DNA2 <- is_revcomp(homodimers$DNA1_sequence[which(is_atcg1 & is_atcg2)], homodimers$DNA2_sequence[which(is_atcg1 & is_atcg2)])

homodimers[which(is_atcg1 & is_atcg2)[which(DNA1_is_rc_of_DNA2==FALSE)], c("PDB ID", paste0("DNA", 1:4, "_sequence"))]

homodimers %>% filter(nakb_DNA3_description!="") %>% pull(`PDB ID`)
#[1] "8R7Z" "4UUV"


metadata <- read_tsv("/projappl/project_2013895/motif_metadata/metadata_final.tsv", col_names=TRUE)

metadata_monomers=metadata %>% filter(experiment %in% c("HT-SELEX", "Methyl-HT-SELEX"))

#better symbols 

TF_info_ensembl= read_delim("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/TF_info_Ensembl115.txt") 
test=left_join(metadata_monomers %>% select(symbol), TF_info_ensembl %>% select(symbol, external_gene_name), by="symbol")

test <- test %>% distinct()   # de-dup exact rows, just in case 719

lst <- test %>%
  filter(!is.na(external_gene_name)) %>%
  group_by(external_gene_name) %>%
  summarise(sym = list(unique(symbol)), .groups = "drop") %>%
  { setNames(.$sym, .$external_gene_name) }

length(lst)


lst[which(sapply(lst, length) >1)]


match_ind=match(metadata_monomers$symbol, test$symbol) 
metadata_monomers$external_gene_name=test$external_gene_name[match_ind]
table(metadata_monomers$external_gene_name==metadata_monomers$symbol)

#139 -> 99
# homodimers=homodimers %>% filter(uniprot1_hgnc_symbol_external_gene_name %in% metadata_monomers$external_gene_name)

TFs=toupper(names(homodimers$uniprot1_hgnc_symbol_external_gene_name %>% table() %>%sort(decreasing=TRUE)))

homodimers=homodimers %>% mutate(motif="") %>% relocate(motif, .after=`DBD type`)
homodimers=homodimers %>% mutate(motif_type="") %>% relocate(motif_type, .after=`motif`)
homodimers=homodimers %>% mutate(seed="") %>% relocate(seed, .after=`motif_type`)
homodimers=homodimers %>% mutate(motif_length=NA) %>% relocate(motif_length, .after=`seed`)
homodimers=homodimers %>% mutate(kmer="") %>% relocate(kmer, .after=`motif_length`)
homodimers=homodimers %>% mutate(orientation="") %>% relocate(orientation, .after=`kmer`)

homodimers=homodimers %>% mutate(motif_match_score=NA) %>% relocate(motif_match_score, .after=`orientation`)
homodimers=homodimers %>% mutate(scaled_motif_match_score=NA) %>% relocate(scaled_motif_match_score, .after=`motif_match_score`)
homodimers=homodimers %>% mutate(motif_match_score_rank=NA) %>% relocate(motif_match_score_rank, .after=`scaled_motif_match_score`)

homodimers=homodimers %>% mutate(SELEX_score=NA) %>% relocate(SELEX_score, .after=`motif_match_score_rank`)
homodimers=homodimers %>% mutate(scaled_SELEX_score=NA) %>% relocate(scaled_SELEX_score, .after=`SELEX_score`)
homodimers=homodimers %>% mutate(SELEX_score_rank=NA) %>% relocate(SELEX_score_rank, .after=`scaled_SELEX_score`)

homodimers=homodimers %>% mutate(jellyfish_score=NA) %>% relocate(jellyfish_score, .after=`SELEX_score_rank`)
homodimers=homodimers %>% mutate(scaled_jellyfish_score=NA) %>% relocate(scaled_jellyfish_score, .after=`jellyfish_score`)
homodimers=homodimers %>% mutate(jellyfish_score_rank=NA) %>% relocate(jellyfish_score_rank, .after=`scaled_jellyfish_score`)

homodimers=homodimers %>% mutate(kmer_alignment_start=NA) %>% relocate(kmer_alignment_start, .after=`orientation`)
homodimers=homodimers %>% mutate(kmer_alignment_end=NA) %>% relocate(kmer_alignment_end, .after=`kmer_alignment_start`)
homodimers=homodimers %>% mutate(structure_sequence_alignment_start=NA) %>% relocate(structure_sequence_alignment_start, .after=`kmer_alignment_end`)
homodimers=homodimers %>% mutate(structure_sequence_alignment_end=NA) %>% relocate(structure_sequence_alignment_end, .after=`structure_sequence_alignment_start`)
homodimers=homodimers %>% mutate(alignment_score=NA) %>% relocate(alignment_score, .after=`structure_sequence_alignment_end`)

library(Biostrings)


source("/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/code/functions.R")


homodimers$motif_exists=FALSE
homodimers$kmers_exists=FALSE
homodimers$match_in_kmer_file=FALSE

# assumes you already have: TFs, homodimers, metadata_monomers, align_with_revcomp()

all_hits <- list()

#RFX1 as uracil in the sequence

for (TF in TFs[1:length(TFs)]) {
  print(TF)
  #TF="NFKB2"
  #TF="NFATC2"
  #TF=TFs[1]
  structure_sequences <- try(homodimers %>%
    filter(toupper(uniprot1_hgnc_symbol_external_gene_name) == TF) %>%
    pull(`DNA1_sequence`) %>%
    unique() %>%
    DNAStringSet(), silent=TRUE)
  
  if(("try-error" %in% class(structure_sequences))){
    next
  }
  
  motifs <- metadata_monomers %>% filter(toupper(external_gene_name) == TF)
  motifs<- motifs %>% filter(experiment=="HT-SELEX")
  message(paste0(TF, ": ", nrow(motifs), " motifs"))
  
  if(nrow(motifs)!=0){
    homodimers$motif_exists[which(toupper(homodimers$uniprot1_hgnc_symbol_external_gene_name)==TF)]=TRUE
  }
  
  hits <- list()  # collect rows for this TF
  
  kmer_file_exists=rep(TRUE, length(motifs$ID))
  names(kmer_file_exists)=motifs$ID
  
  for (motif in motifs$ID) {
    #motif=motifs$ID[1]
    message(motif)
    
    
    kmers_scores=try(read_delim(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",motif,"_combined_counts.tsv"), col_names = TRUE),
                     silent=TRUE)
    
    if(("try-error" %in% class(kmers_scores))){
      #kmers_scores_orig=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/",motif,"_score_rank.tsv"), col_names = FALSE)
      kmers_scores=try(read_delim(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",motif,"_jellyfish_counts.tsv"), col_names = TRUE),
                       silent=TRUE)
      
    }
    
    if(("try-error" %in% class(kmers_scores))){
      #kmers_scores=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/",motif,".tsv"), col_names = FALSE)
      #kmers_scores_orig=read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/",motif,".tsv"), col_names = FALSE)
      
      kmers_scores=try(read_delim(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",motif,".tsv"), col_names = FALSE),
                       silent=TRUE)
      
    }
    if(("try-error" %in% class(kmers_scores))){
      kmer_file_exists[motif]=FALSE
      next
    }
    
    if(ncol(kmers_scores)==4){
    
      names(kmers_scores)=c("kmers", "motif_match_score","scaled_motif_match_score", "motif_match_score_rank")
      #
    }else if(ncol(kmers_scores)==29){ #combined
      names(kmers_scores)[1:4]=c("kmers", "motif_match_score","scaled_motif_match_score", "motif_match_score_rank")
    }else if(ncol(kmers_scores)==18){ #SELEX only
      names(kmers_scores)[1:4]=c("kmers", "motif_match_score","scaled_motif_match_score", "motif_match_score_rank")
    }else if(ncol(kmers_scores)==15){ #jellyfish only
      names(kmers_scores)[1:4]=c("kmers", "motif_match_score","scaled_motif_match_score", "motif_match_score_rank")
    }else{
      print("wrong number of columns")
    }
    
    
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
      
      make_kmer_cols <- function(kmers_scores, idx, missing_as_empty = FALSE) {
        # map: source name in kmers_scores  ->  output column name
        src_to_new <- c(
          "motif_match_score"                = "motif_match_score",
          "scaled_motif_match_score"         = "scaled_motif_match_score",
          "motif_match_score_rank"           = "motif_match_score_rank",
          "Corrected count"                  = "SELEX_score",
          "Corrected count scaled"           = "scaled_SELEX_score",
          "Corrected_count_rank"             = "SELEX_score_rank",
          "jellyfish Corrected count"        = "jellyfish_score",
          "jellyfish Corrected count scaled" = "scaled_jellyfish_score",
          "jellyfish Corrected_count_rank"   = "jellyfish_score_rank"
          
         
          
          
        )
        
        # preferred defaults (typed NAs)
        defaults_na <- list(
          motif_match_score         = NA_real_,
          scaled_motif_match_score  = NA_real_,
          motif_match_score_rank    = NA_integer_,
          SELEX_score               = NA_real_,
          scaled_SELEX_score        = NA_real_,
          SELEX_score_rank          = NA_integer_,
          jellyfish_score           = NA_real_,
          scaled_jellyfish_score    = NA_real_,
          jellyfish_score_rank      = NA_integer_
        )
        # or: fill with empty strings instead of NA (will make those columns character)
        defaults <- if (missing_as_empty) lapply(defaults_na, function(.) "") else defaults_na
        
        # start with all defaults present
        vals <- defaults
        
        # overwrite from kmers_scores when the source column exists
        present_src <- intersect(names(src_to_new), names(kmers_scores))
        for (src in present_src) {
          new <- src_to_new[[src]]
          v   <- kmers_scores[[src]][idx]
          if (length(v) == 0) v <- defaults[[new]]
          vals[[new]] <- v
        }
        vals
      }
      
    
      
      
      
      # Collect a tidy row
      hits[[length(hits) + 1L]] <- tibble(
        TF = TF,
        motif_id = motif,
        motif_length = metadata_monomers %>% filter(ID == motif) %>% pull(length),
        seed        = metadata_monomers %>% filter(ID == motif) %>% pull(seed),
        motif_type  = metadata_monomers %>% filter(ID == motif) %>% pull(type),
        structure_sequence = as.character(subject),
        kmer              = as.character(kmers[idx]),
        orientation       = orientation,
        alignment_score   = score(aln),
        kmer_alignment_start = start(aln@pattern@range),
        kmer_alignment_end   = end(aln@pattern@range),
        structure_sequence_alignment_start = start(aln@subject@range),
        structure_sequence_alignment_end   = end(aln@subject@range),
        !!!make_kmer_cols(kmers_scores, idx, missing_as_empty = FALSE)  # TRUE -> fill "" instead of NA
      )
        
        
        # tibble(
        # TF = TF,
        # motif_id = motif,
        # motif_length=metadata_monomers %>% filter(ID==motif) %>% pull(length),
        # seed=metadata_monomers %>% filter(ID==motif) %>% pull(seed),
        # motif_type=metadata_monomers %>% filter(ID==motif) %>% pull(type),
        # structure_sequence = as.character(subject),
        # kmer = as.character(kmers[idx]),
        # orientation = orientation,
        # alignment_score = score(aln),
        # kmer_alignment_start=start(aln@pattern@range), #1
        # kmer_alignment_end=end(aln@pattern@range), #15
        # structure_sequence_alignment_start= start(aln@subject@range), #4
        # structure_sequence_alignment_end=end(aln@subject@range),
        # motif_match_score = kmers_scores$motif_match_score[idx],
        # motif_match_score_scaled = kmers_scores$motif_match_score_scaled[idx],
        # motif_match_score_rank = kmers_scores$motif_match_rank[idx],
        # SELEX_score=kmers_scores$`Corrected count`[idx],
        # scaled_SELEX_score=kmers_scores$`Corrected count scaled`[idx],
        # SELEX_score_rank=kmers_scores$Corrected_count_rank[idx],
        # jellyfish_score=kmers_scores$`jellyfish Corrected count`[idx],
        # scaled_jellyfish_score=kmers_scores$`jellyfish Corrected count scaled`[idx],
        # jellyfish_score_rank=kmers_scores$`jellyfish Corrected_count_rank`[idx]
        # )  #18
        
    } #sequences
  } #motifs
  
  if(length(grep("TRUE", names(table(kmer_file_exists))))>0){
    homodimers$kmers_exists[which(toupper(homodimers$uniprot1_hgnc_symbol_external_gene_name)==TF)]=TRUE
  }
  
  all_hits[[TF]] <- if (length(hits)) bind_rows(hits) else tibble()

  if(nrow( all_hits[[TF]])==0){
    homodimers$match_in_kmer_file[which(toupper(homodimers$uniprot1_hgnc_symbol_external_gene_name)==TF)]=FALSE
  }

  if(nrow(all_hits[[TF]])>0){
    df=all_hits[[TF]]  %>% group_by(structure_sequence) %>% nest()
    
    
    for(j in 1:length(df$structure_sequence)){
      #j=1
      seq=df$structure_sequence[j]
      max_ind=which(df$data[[j]]$alignment_score==max(df$data[[j]]$alignment_score))
      
      df$data[[j]]=df$data[[j]][max_ind,]
      if(unique(!is.na(df$data[[j]]$SELEX_score_rank))){
        min_ind=which(df$data[[j]]$SELEX_score_rank==min(df$data[[j]]$SELEX_score_rank))
      }else{
        min_ind=which(df$data[[j]]$motif_match_score_rank==min(df$data[[j]]$motif_match_score_rank))
      }
      
      df$data[[j]]=df$data[[j]][min_ind,]
      
      if(nrow(df$data[[j]])>1){
        print("multiple matches")
      }
     
      ind=which(toupper(homodimers$uniprot1_hgnc_symbol_external_gene_name)==TF & homodimers$`DNA1_sequence`==seq)
      print(ind)
      homodimers$match_in_kmer_file[ind]=TRUE
      
      homodimers$motif[ind]=df$data[[j]]$motif_id
      homodimers$motif_type[ind]=df$data[[j]]$motif_type
      homodimers$seed[ind]=df$data[[j]]$seed
      homodimers$motif_length[ind]=df$data[[j]]$motif_length
      homodimers$kmer[ind]=df$data[[j]]$kmer
      homodimers$orientation[ind]=df$data[[j]]$orientation
  
      
      homodimers$kmer_alignment_start[ind]=df$data[[j]]$kmer_alignment_start
      homodimers$kmer_alignment_end[ind]=df$data[[j]]$kmer_alignment_end
      homodimers$structure_sequence_alignment_start[ind]= df$data[[j]]$structure_sequence_alignment_start
      homodimers$structure_sequence_alignment_end[ind]=df$data[[j]]$structure_sequence_alignment_end
      homodimers$alignment_score[ind]=df$data[[j]]$alignment_score
      
      
      
      homodimers$motif_match_score[ind]=df$data[[j]]$motif_match_score
      homodimers$scaled_motif_match_score[ind]=df$data[[j]]$scaled_motif_match_score
      homodimers$motif_match_score_rank[ind]=df$data[[j]]$motif_match_score_rank
      
      
      homodimers$SELEX_score[ind]=df$data[[j]]$SELEX_score
      homodimers$scaled_SELEX_score[ind]=df$data[[j]]$scaled_SELEX_score
      homodimers$SELEX_score_rank[ind]=df$data[[j]]$SELEX_score_rank
      
      homodimers$jellyfish_score[ind]=df$data[[j]]$jellyfish_score
      homodimers$scaled_jellyfish_score[ind]=df$data[[j]]$scaled_jellyfish_score
      homodimers$jellyfish_score_rank[ind]=df$data[[j]]$jellyfish_score_rank
      
    } #over sequences, find the best motif
  
  } #hits were found

} #over TFs


homodimers=homodimers %>% relocate(uniprot1_hgnc_symbol_external_gene_name, .after=TFs)

homodimers=homodimers %>% relocate(match_in_kmer_file, .after=`DBD type`)
homodimers=homodimers %>% relocate(kmers_exists, .after=`DBD type`)
homodimers=homodimers %>% relocate(motif_exists, .after=`DBD type`)

raw="
 [1] \"TFs\"                                     \"uniprot1_hgnc_symbol_external_gene_name\" \"TFs_alternative\"                         \"TF_note\"                                
 [5] \"PDB ID\"                                  \"DBD type\"                        
 [4] \"motif_exists\"                            \"kmers_exists\"                           
 [9] \"match_in_kmer_file\"       \"motif\"                                   \"motif_type\"                             
 [9] \"seed\"                                    \"motif_length\"                            \"kmer\"                                    \"orientation\"                            
[13] \"kmer_alignment_start\"                    \"kmer_alignment_end\"                      \"structure_sequence_alignment_start\"      \"structure_sequence_alignment_end\"       
[17] \"alignment_score\"                         \"motif_match_score\"                       \"scaled_motif_match_score\"                \"motif_match_score_rank\"                 
[21] \"SELEX_score\"                             \"scaled_SELEX_score\"                      \"SELEX_score_rank\"                        \"jellyfish_score\"                        
[25] \"scaled_jellyfish_score\"                  \"jellyfish_score_rank\"                    \"Resolution (‚Ñ´)\"                          \"DNA sequence\"                           
[29] \"notes\"                                   \"TF-TF contact area\"                      \"Colour legend\"                           \"Type\"                                   
[33] \"DNA_count\"                               \"Protein_count\"                           \"DNA1_sequence\"                           \"DNA2_sequence\"                          
[37] \"DNA3_sequence\"                           \"DNA4_sequence\"                           \"Protein1_sequence\"                       \"Protein2_sequence\"                      
[41] \"Protein3_sequence\"                       \"Protein4_sequence\"                       \"Protein1_molecule\"                       \"Protein2_molecule\"                      
[45] \"Protein3_molecule\"                       \"Protein4_molecule\"                       \"DNA1_organism\"                           \"DNA2_organism\"                          
[49] \"DNA3_organism\"                           \"DNA4_organism\"                           \"Protein1_organism\"                       \"Protein2_organism\"                      
[53] \"Protein3_organism\"                       \"Protein4_organism\"                       \"ref1\"                                    \"ref2\"                                   
[57] \"ref3\"                                    \"ref4\"                                    \"ref5\"                                    \"ref6\"                                   
[61] \"struct_title\"                            \"struct_pdbx_keywords\"                    \"struct_text_keywords\"                    \"NAKBpolyclass\"                          
[65] \"NAKBmethod\"                              \"NAKBna\"                                  \"NAKBprot\"                                \"uniprot1\"                               
[69] \"uniprot2\"                                \"uniprot3\"                                \"uniprot4\"                                \"nakb_Protein1_uniprot\"                  
[73] \"nakb_Protein2_uniprot\"                   \"nakb_Protein3_uniprot\"                   \"nakb_Protein4_uniprot\"                   \"uniprot2_hgnc_symbol_external_gene_name\"
[77] \"uniprot3_hgnc_symbol_external_gene_name\" \"uniprot4_hgnc_symbol_external_gene_name\" "

cols <- str_match_all(raw, "\"([^\"]+)\"")[[1]][,2]
cols <- cols[nzchar(cols)]
cols <- cols[!duplicated(cols)]  # keep order, drop dups

test=homodimers[,cols]


write_delim(test,file="/projappl/project_2013895/SELEX/PDB_structures_to_motifs_and_kmers/Data/homodimers_better.tsv", delim="\t", col_names=TRUE)  
    
  
  
  
  