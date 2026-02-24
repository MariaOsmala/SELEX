library("readr")
library("dplyr")
library("tidyverse")
rm(list=ls())
metadata=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv", delim="\t")

which(metadata$ID=="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3")

count_kmers_with_N=read_delim("/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3_hamming1_signal.txt",
                              col_names = FALSE)

names(count_kmers_with_N)=c("pattern", "deg_count")

#jellyfish_counts=read_delim("/scratch/project_2013895/SELEX/jellyfish/filtered_counts/ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3_signal_counts.txt",
#                           col_names=FALSE)

jellyfish_counts=read_delim("/scratch/project_2013895/SELEX/jellyfish/counts/ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3_signal_counts.txt",
                           col_names=FALSE)

names(jellyfish_counts)=c("kmer", "count")

#convert N-pattern to regex
pattern_to_regex <- function(pat) {
  regex_body <- gsub("N", ".", pat)
  paste0("^", regex_body, "$")
}

#We’ll convert each pattern, e.g. CACGTGNNNGGATTAN → ^CACGTG...GGATTA.$
# Reverse complement, vectorized and safe for character(0)
revcomp <- function(s) {
  if (length(s) == 0) return(s)  # character(0) → character(0)
  
  sapply(s, function(one) {
    chars <- strsplit(one, "")[[1]]
    comp  <- chartr("ACGTacgt", "TGCAtgca", paste(chars, collapse = ""))
    rev_chars <- rev(strsplit(comp, "")[[1]])
    paste(rev_chars, collapse = "")
  }, USE.NAMES = FALSE)
}




jellyfish_counts <- jellyfish_counts%>%
  mutate(rev_comp = revcomp(kmer))



#3. For each degenerate k-mer, sum counts of matching exact k-mers
#This uses rowwise() so each pattern is processed independently:

count_kmers_with_N[21,]
#CACGTGNNNAGATTAN       238
pattern=count_kmers_with_N$pattern[21]
regex = pattern_to_regex(pattern)
matches=grepl(regex, jellyfish_counts$kmer)
match_seqs=jellyfish_counts$kmer[matches] #192
sum_exact_counts = sum(jellyfish_counts$count[matches])

matches=grepl(revcomp(regex), jellyfish_counts$kmer)
match_seqs=jellyfish_counts$kmer[matches] #192
sum_exact_counts_revcomp = sum(jellyfish_counts$count[matches])

sum_exact_counts+sum_exact_counts_revcomp

result <- count_kmers_with_N %>%
  rowwise() %>%
  mutate(
    #pattern=count_kmers_with_N$pattern[21]
    regex = pattern_to_regex(pattern),
    
    # forward matches
    matches_fwd = list(grepl(regex, jellyfish_counts$kmer)),
    sum_exact_counts = sum(jellyfish_counts$count[unlist(matches_fwd)]),
    
    # reverse-complement matches
    matches_rev = list(grepl(regex, jellyfish_counts$rev_comp)),
    sum_exact_counts_revcomp = sum(jellyfish_counts$count[unlist(matches_rev)]),
    
    total_sum = sum_exact_counts + sum_exact_counts_revcomp
  ) %>%
  ungroup() %>%
  select(-matches_fwd, -matches_rev)  # drop the logical list-columns if you don't need them

(result$deg_count==result$total_sum) %>% table()

#TRUE
#70
test=result[which(!(result$deg_count==result$total_sum)),] #empty

