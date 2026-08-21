library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("tidyverse")
library(R.utils)
library(stringr)
library(tidyr)
library("ggpubr")


rm(list=ls())

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635
metadata=metadata[order(metadata$length),]

correlation_data=data.frame(ID=metadata$ID)
correlation_data$pearson_estimate=NA
correlation_data$spearman_estimate=NA
correlation_data$kendall_estimate=NA

correlation_data$pearson_pvalue=NA
correlation_data$spearman_pvalue=NA
correlation_data$kendall_pvalue=NA

figure_path="/scratch/project_2013895/SELEX/Figures/correlation_between_motif_match_scores_and_kmer_counts/"

all_data=data.frame()

for(i in 1:nrow(metadata)){
  print(i)
  #i=1
  if(file.exists(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",metadata$ID[i],"_combined_counts.tsv"))){
    file=paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",metadata$ID[i],"_combined_counts.tsv")
  }else if(file.exists(paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",metadata$ID[i],"_jellyfish_counts.tsv"))){
    file=paste0("/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/",metadata$ID[i],"_combined_counts.tsv")
  }else{
    file=NULL
    next
  }
     
   data=read_delim(file=file, delim="\t")

   if(unique(is.na(data$`jellyfish Corrected count scaled`))){
     next
   }
   
   all_data=rbind(all_data, data %>% select(motif_match_score_scaled, `jellyfish Corrected count scaled`))
     
   #cor(data$motif_match_score_scaled, data$`jellyfish Corrected count scaled`, method="pearson")
   #cor(data$motif_match_score_scaled, data$`jellyfish Corrected count scaled`, method="spearman")
   
   ct_pearson=cor.test(data$motif_match_score_scaled, data$`jellyfish Corrected count scaled`, method="pearson", alternative="two.sided")
   ct_spearman=cor.test(data$motif_match_score_scaled, data$`jellyfish Corrected count scaled`, method="spearman", alternative="two.sided")
   ct_kendall=cor.test(data$motif_match_score_scaled, data$`jellyfish Corrected count scaled`, method="kendall", alternative="two.sided")
   
   correlation_data$pearson_estimate[i]=as.numeric(ct_pearson$estimate)
   correlation_data$spearman_estimate[i]=as.numeric(ct_spearman$estimate)
   correlation_data$kendall_estimate[i]=as.numeric(ct_kendall$estimate)
   
   correlation_data$pearson_pvalue[i]=as.numeric(ct_pearson$p.value)
   correlation_data$spearman_pvalue[i]=as.numeric(ct_spearman$p.value)
   correlation_data$kendall_pvalue[i]=as.numeric(ct_kendall$p.value)
   
   p=ggscatter(data, x = "motif_match_score_scaled", y = "jellyfish Corrected count scaled", 
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "spearman",
             xlab = "Scaled motif match score", ylab = "Scaled k-mer count score")
   

   ggsave(
     filename = paste0(figure_path, metadata$ID[i], ".pdf"),
     plot = p,
     device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
     width = 7, height = 5, units = "in"  # set page size
   )          
   
}

saveRDS(correlation_data, "/scratch/project_2013895/SELEX/RData/correlation_between_motif_match_scores_and_count_scores_correlations.RDS")

#correlation_data=readRDS("/scratch/project_2013895/SELEX/RData/correlation_between_motif_match_scores_and_count_scores_correlations.RDS")

saveRDS(all_data, "/scratch/project_2013895/SELEX/RData/correlation_between_motif_match_scores_and_count_scores_all_data.RDS")

#all_data=readRDS("/scratch/project_2013895/SELEX/RData/correlation_between_motif_match_scores_and_count_scores_all_data.RDS")

#1. Fast & simple: plot a random subset
#You almost never need all 18M points to see the pattern.

set.seed(1)
n_sample <- 500000  # e.g. 100k points, You can tune n_sample (e.g. 50k–500k) to balance speed vs detail.
idx <- sample.int(nrow(all_data), n_sample)

plot(all_data$motif_match_score_scaled[idx], all_data$`jellyfish Corrected count scaled`[idx],
     pch = ".", cex = 1,
     xlab = "x", ylab = "y")

#Use smoothScatter (base R, density-based)
#This draws a density image instead of every point.
library(viridis)
library(fields)
png(paste0(figure_path, "all_data_density_plot.png"), width = 2000, height = 2000, res = 300)
with(all_data%>% filter((motif_match_score_scaled>=-1 )& ( `jellyfish Corrected count scaled`>=-1)), smoothScatter(motif_match_score_scaled, `jellyfish Corrected count scaled`,
                             xlim=c(-1,1), ylim=c(-1,1),
                             nbin = 512,              # resolution of the grid, the number of equally spaced grid points for the density estimation;
                             bandwidth = 0.01,        # smoothing bandwidth, the bandwidth to be used in each coordinate direction.
                             #gridsize=128,
                             nrpoints = 0,            # how many raw points to overplot, number of points to be superimposed on the density image. The first nrpoints points from those areas of lowest regional densities will be plotted. Adding points to the plot allows for the identification of outliers. If all points are to be plotted, choose nrpoints = Inf.
                             colramp = colorRampPalette(c("white", "blue", "red")) ,
                             #colramp       = viridis, 
     #                        transformation = function(z) log1p(z),  # transform counts
                       xlab = "Motif match scaled score",
                       ylab = "k-mer count (jellyfish) scaled score"))

dev.off()


#Hexbin / 2D histogram (aggregation)
#Using hexbin package
#install.packages("hexbin")  # if needed
library(hexbin)

# This gives a hexagon-aggregated view where color = point density.
png(paste0(figure_path, "all_data_hexbin.png"), width = 2000, height = 2000, res = 300)
with(all_data, {
  h <- hexbin(motif_match_score_scaled, `jellyfish Corrected count scaled`, xbins = 100,
              xlab = "Motif match scaled score",
              ylab = "k-mer count (jellyfish) scaled score" )  # adjust xbins for resolution
  plot(h)
})
dev.off()

with(all_data %>% filter((motif_match_score_scaled>=-5 )& ( `jellyfish Corrected count scaled`>=-1)), {
  h <- hexbin(motif_match_score_scaled, `jellyfish Corrected count scaled`, xbins = 500,
              xlab = "Motif match scaled score",
              ylab = "k-mer count (jellyfish) scaled score" )  # adjust xbins for resolution
  plot(h)
})


ggplot(all_data%>% filter((motif_match_score_scaled>=-5 )& ( `jellyfish Corrected count scaled`>=-1)), aes(motif_match_score_scaled, `jellyfish Corrected count scaled`)) +
  # geom_bin2d() +  # or 
  stat_bin_2d() +
  labs( x= "Motif match scaled score",
        y = "k-mer count (jellyfish) scaled score")

ggplot(all_data, aes(motif_match_score_scaled, `jellyfish Corrected count scaled`)) +
  geom_binhex() +  # or stat_bin_2d()
  labs(x = "x", y = "y")

p=ggscatter(all_data%>% filter((motif_match_score_scaled>=-5 )& ( `jellyfish Corrected count scaled`>=-1)), x = "motif_match_score_scaled", y = "jellyfish Corrected count scaled", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Scaled motif match score", ylab = "Scaled k-mer count score")


ggsave(
  filename = paste0(figure_path, "all_data.pdf"),
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)          




