library("readr")
library("dplyr")
library("ggpubr")

rm(list=ls())

# bHLH_Homeodomain pair ARNTL_PITX
motif="ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3"

#Protein sequences of this motif

proteins=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv")

data=proteins %>% filter(ID==motif)
data$`protein sequence 3`=data$`protein sequence 2`
data$`protein sequence 2`=data$`protein sequence 1`

data=data[,-c(12,13,14,15)]
data$type=NULL
data$clone=NULL
write_delim(data, paste0("/projappl/project_2013895/SELEX/metadata.tsv"), delim="\t")

#all realisations and Hamming1 neighbourhood but substitutions only at non-N positions
all_realisations=read_delim(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",motif,"_combined_counts.tsv")) 

names(all_realisations)

correlation_data=data.frame(ID=motif)
correlation_data$pearson_estimate=NA
correlation_data$spearman_estimate=NA
correlation_data$kendall_estimate=NA

correlation_data$pearson_pvalue=NA
correlation_data$spearman_pvalue=NA
correlation_data$kendall_pvalue=NA


ct_pearson=cor.test(all_realisations$motif_match_score_scaled, all_realisations$`jellyfish Corrected count scaled`, method="pearson", alternative="two.sided")
ct_spearman=cor.test(all_realisations$motif_match_score_scaled, all_realisations$`jellyfish Corrected count scaled`, method="spearman", alternative="two.sided")
ct_kendall=cor.test(all_realisations$motif_match_score_scaled, all_realisations$`jellyfish Corrected count scaled`, method="kendall", alternative="two.sided")
i=1
correlation_data$pearson_estimate[i]=as.numeric(ct_pearson$estimate)
correlation_data$spearman_estimate[i]=as.numeric(ct_spearman$estimate)
correlation_data$kendall_estimate[i]=as.numeric(ct_kendall$estimate)

correlation_data$pearson_pvalue[i]=as.numeric(ct_pearson$p.value)
correlation_data$spearman_pvalue[i]=as.numeric(ct_spearman$p.value)
correlation_data$kendall_pvalue[i]=as.numeric(ct_kendall$p.value)

p=ggscatter(all_realisations, x = "motif_match_score_scaled", y = "jellyfish Corrected count scaled", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Scaled motif match score", ylab = "Scaled k-mer count score")


kmers_with_N=read_delim(paste0("/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/", motif, ".tsv"))

correlation_data=data.frame(ID=motif)
correlation_data$pearson_estimate=NA
correlation_data$spearman_estimate=NA
correlation_data$kendall_estimate=NA

correlation_data$pearson_pvalue=NA
correlation_data$spearman_pvalue=NA
correlation_data$kendall_pvalue=NA


ct_pearson=cor.test(kmers_with_N$motif_match_score_scaled, kmers_with_N$Corrected_count_scaled, method="pearson", alternative="two.sided")
ct_spearman=cor.test(kmers_with_N$motif_match_score_scaled, kmers_with_N$Corrected_count_scaled, method="spearman", alternative="two.sided")
ct_kendall=cor.test(kmers_with_N$motif_match_score_scaled, kmers_with_N$Corrected_count_scaled, method="kendall", alternative="two.sided")
i=1
correlation_data$pearson_estimate[i]=as.numeric(ct_pearson$estimate)
correlation_data$spearman_estimate[i]=as.numeric(ct_spearman$estimate)
correlation_data$kendall_estimate[i]=as.numeric(ct_kendall$estimate)

correlation_data$pearson_pvalue[i]=as.numeric(ct_pearson$p.value)
correlation_data$spearman_pvalue[i]=as.numeric(ct_spearman$p.value)
correlation_data$kendall_pvalue[i]=as.numeric(ct_kendall$p.value)

p=ggscatter(kmers_with_N, x = "motif_match_score_scaled", y = "Corrected_count_scaled", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Scaled motif match score", ylab = "Scaled k-mer count score")





all_realisations =all_realisations %>% select(c("Kmer","Corrected count scaled" ,"Corrected_count_rank"))

all_realisations$`Corrected count scaled`=as.double(all_realisations$`Corrected count scaled`)

all_realisations=all_realisations[order(all_realisations$`Corrected count scaled`, decreasing=TRUE),]

all_realisations$Corrected_count_rank=1:nrow(all_realisations)

all_realisations <- all_realisations %>% rename(Corrected_count_scaled = `Corrected count scaled`)

data.table::setnames(all_realisations, old = "Corrected count scaled", new = "Corrected_count_scaled ")

write_delim(all_realisations, paste0("/projappl/project_2013895/SELEX/",motif,"_without_Ns.tsv"), delim="\t")

#realisations with N keeped, Hamming1 and Hamming2 neighbourhood substutions only at non-N positions
N_keeped=read_delim("/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/ARNTL_PITX1_CAP-SELEX_TCTTTC40NTTG_AAC_CACGTGNNNRGATTAN_1_3.tsv")

#Which correspond to what


names(N_keeped)=c("pattern", "deg_count", "deg_count_rank")


names(all_realisations)=c("kmer", "count")

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




all_realisations <- all_realisations%>%
  mutate(rev_comp = revcomp(kmer))



#3. For each degenerate k-mer, sum counts of matching exact k-mers
#This uses rowwise() so each pattern is processed independently:

count_kmers_with_N[21,]
#CACGTGNNNAGATTAN       238
pattern=count_kmers_with_N$pattern[21]
regex = pattern_to_regex(pattern)
matches=grepl(regex, all_realisations$kmer)
match_seqs=all_realisations$kmer[matches] #192
sum_exact_counts = sum(all_realisations$count[matches])

matches=grepl(revcomp(regex), all_realisations$kmer)
match_seqs=all_realisations$kmer[matches] #192
sum_exact_counts_revcomp = sum(all_realisations$count[matches])

sum_exact_counts+sum_exact_counts_revcomp

result <- count_kmers_with_N %>%
  rowwise() %>%
  mutate(
    #pattern=count_kmers_with_N$pattern[21]
    regex = pattern_to_regex(pattern),
    
    # forward matches
    matches_fwd = list(grepl(regex, all_realisations$kmer)),
    sum_exact_counts = sum(all_realisations$count[unlist(matches_fwd)]),
    
    # reverse-complement matches
    matches_rev = list(grepl(regex, all_realisations$rev_comp)),
    sum_exact_counts_revcomp = sum(all_realisations$count[unlist(matches_rev)]),
    
    total_sum = sum_exact_counts + sum_exact_counts_revcomp
  ) %>%
  ungroup() %>%
  select(-matches_fwd, -matches_rev)  # drop the logical list-columns if you don't need them

(result$deg_count==result$total_sum) %>% table()

#TRUE
#70
test=result[which(!(result$deg_count==result$total_sum)),] #empty





