library("ggplot2")
library(readr)
library(tidyverse)
library(dplyr)


successfull=data.frame(ID=gsub("_score_rank.tsv","",
                 dir("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore")
))


info=data.frame(ID=rep(NA,nrow(successfull)), 
                kmer_nro=rep(NA,nrow(successfull)), 
                min_score=rep(NA,nrow(successfull)), 
                max_score=rep(NA,nrow(successfull)), 
                min_scaled_score=rep(NA,nrow(successfull)), 
                max_scaled_score=rep(NA,nrow(successfull)))

hist_path="/scratch/project_2013895/SELEX/Hamming12_score_histograms/"

for(i in 1:nrow(successfull)){
  #i=1
  print(i)
  info$ID[i]=successfull$ID[i]
  data=try(read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/", 
                             successfull$ID[i],"_score_rank.tsv"), 
                      delim="\t", col_names=FALSE), silent=TRUE)
  if( unique(class(data)=="try-error" ) ){
    next
  }else if( ncol(data)!=3){
    next
  }
  names(data)=c("kmer", "score", "rank")
  
 
  p <- ggplot(data, aes(x = score)) +
    geom_histogram() +
    scale_x_continuous() +
    labs(
      title = successfull$ID[i],
      x = "score",
      y = "count"
    ) +
    theme_minimal() +theme(plot.title = element_text(size = 10) )
  
  #print(p)
  ggsave(paste0(hist_path, "scores/",successfull$ID[i],".png"), 
                plot = p, width = 6, height = 4, dpi = 300, units = "in", bg = "white")
  
 
  info$min_score[i]=min(data$score)
  info$max_score[i]=max(data$score)
  info$kmer_nro[i]=nrow(data)
  
  data=try(read_delim(paste0("/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/", 
                             successfull$ID[i],"_score_rank.tsv"),delim="\t",col_names=FALSE), silent=TRUE)
  
  if( unique( class(data)=="try-error" ) ){
    next
  }else if( ncol(data)!=3 ){
    next
  }
  
  names(data)=c("kmer", "score", "rank")
  
  p <- ggplot(data, aes(x = score)) +
    geom_histogram() +
    scale_x_continuous() +
    labs(
      title = successfull$ID[i],
      x = "score",
      y = "count"
    ) +
    theme_minimal() +theme(plot.title = element_text(size = 10) )
  
  #print(p)
  ggsave(paste0(hist_path, "scaled_scores/",successfull$ID[i],".png"), 
         plot = p, width = 6, height = 4, dpi = 300, units = "in", bg = "white")
  
  info$min_scaled_score[i]=min(data$score)
  info$max_scaled_score[i]=max(data$score)

}

write.table(info, 
            "/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming12_score_analysis.tsv",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)


stop()

info_hamming12=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming12_score_analysis.tsv")

info_hamming1=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/score_analysis.tsv")

info_hamming1_updated=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming1_score_analysis.tsv")

names(info_hamming1_updated)=c("ID", paste0("Hamming1_updated_", names(info_hamming1_updated)[-1]))

names(info_hamming1)=c("ID", paste0("Hamming1_", names(info_hamming1)[-1]))
names(info_hamming12)=c("ID", paste0("Hamming12_", names(info_hamming12)[-1]))

test=info_hamming12 %>% left_join(info_hamming1, by = "ID")

test=test %>% left_join(info_hamming1_updated, by="ID")

test=test[,c("ID", "Hamming1_kmer_nro", "Hamming1_updated_kmer_nro","Hamming12_kmer_nro")]


ggplot(test, aes(x = Hamming1_updated_kmer_nro)) +
  geom_histogram(bins=100) + 
  scale_x_continuous() + xlim(c(0,1000))+
  labs(
    title = "Hamming1",
    x = "k-mer number",
    y = "count"
  ) +
  theme_minimal() +theme(plot.title = element_text(size = 10) )


ggplot(test, aes(x = Hamming12_kmer_nro)) +
  geom_histogram(bins=100) + 
  scale_x_continuous() + xlim(c(0,1000))+
  labs(
    title = "Hamming1",
    x = "k-mer number",
    y = "count"
  ) +
  theme_minimal() +theme(plot.title = element_text(size = 10) )



metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")


test=test %>% left_join(metadata %>% select(ID, length, IC)

ggplot(test, aes(x = length, y = Hamming1_updated_kmer_nro, color=IC/length)) +
  scale_color_gradient(low = "blue", high = "red")+
  geom_jitter(width = 0.2, height = 0.2)+
  scale_y_log10() +
  labs(title = "Hamming1", x = "seed length", y = "k-mer number") +
  theme_minimal()


ggplot(test, aes(x = length, y = Hamming12_kmer_nro, color=IC/length)) +
  scale_color_gradient(low = "blue", high = "red")+
  geom_jitter(width = 0.2, height = 0.2)+
  scale_y_log10() +
  labs(title = "Hamming12", x = "seed length", y = "k-mer number") +
  theme_minimal()




