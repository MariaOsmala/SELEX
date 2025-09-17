library("ggplot2")
library(readr)
library(tidyverse)
library(dplyr)



successfull=data.frame(ID=gsub(".tsv","",dir("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore")
)
)


info=data.frame(ID=rep(NA,nrow(successfull)), 
                kmer_nro=rep(NA,nrow(successfull)), 
                min_score=rep(NA,nrow(successfull)), 
                max_score=rep(NA,nrow(successfull)), 
                min_scaled_score=rep(NA,nrow(successfull)), 
                max_scaled_score=rep(NA,nrow(successfull)))

hist_path="/scratch/project_2013895/SELEX/Hamming1_score_histograms/"

for(i in 1:nrow(successfull)){
  #i=1
  print(i)
  info$ID[i]=successfull$ID[i]
  data=try(read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/", successfull$ID[i],".tsv"), delim="\t", col_names=FALSE), silent=TRUE)
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
  
  data=try(read_delim(paste0("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/", successfull$ID[i],".tsv"),delim="\t", col_names=FALSE), silent=TRUE)
  
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
            "/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming1_score_analysis.tsv",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)




