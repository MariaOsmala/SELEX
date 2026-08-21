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


#metadata=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv", delim="\t")

#For how may motifs, SELEX or background data missing

table(is.na(metadata$CSC_SELEX_filename) | is.na(metadata$CSC_SELEX_background_filename))
#TRUE 62

#For how many lambda is missing

correlation_data=data.frame(ID=metadata$ID)
correlation_data$pearson_estimate=NA
correlation_data$spearman_estimate=NA
correlation_data$kendall_estimate=NA

correlation_data$pearson_pvalue=NA
correlation_data$spearman_pvalue=NA
correlation_data$kendall_pvalue=NA

figure_path="/scratch/project_2013895/SELEX/Figures/correlation_between_motif_match_scores_and_kmer_counts/"

results_path="/scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds/"

all_data=data.frame()

for(i in 1:nrow(metadata)){
  print(i)
  #i=1
  #i=708
  if(file.exists(paste0(results_path,metadata$ID[i],".tsv"))){
    file=paste0(results_path,metadata$ID[i],".tsv")
  }else{
    file=NULL
    next
  }
     
   data=read_delim(file=file, delim="\t")
   if( names(table(is.na(data$Corrected_count_scaled)))=="TRUE"){
     next
   }
   
   data=data %>% mutate(ID=metadata$ID[i]) %>% relocate(ID, .before=Kmer)  
   all_data=rbind(all_data, data)
   
   #measures linear correlation between two sets of data.
   #the ratio between the covariance of two variables and the product of their standard deviations; 
   #thus, it is essentially a normalized measurement of the covariance, such that the result always has a value between −1 and 1.
   ct_pearson=cor.test(data$motif_match_score_scaled, data$Corrected_count_scaled, method="pearson", alternative="two.sided")
   #How strongly two sets of ranks are correlated. Whether k-mers who are high ranking in SELEX counts are also high ranking in motif matching
   #A nonparametric measure of rank correlation (statistical dependence between the rankings of two variables). 
   #It assesses how well the relationship between two variables can be described using a monotonic function.
   ct_spearman=cor.test(data$motif_match_score_scaled, data$Corrected_count_scaled, method="spearman", alternative="two.sided")
   #a statistic used to measure the ordinal association between two measured quantities.
   ct_kendall=cor.test(data$motif_match_score_scaled, data$Corrected_count_scaled, method="kendall", alternative="two.sided")
   
   correlation_data$pearson_estimate[i]=as.numeric(ct_pearson$estimate)
   correlation_data$spearman_estimate[i]=as.numeric(ct_spearman$estimate)
   correlation_data$kendall_estimate[i]=as.numeric(ct_kendall$estimate)
   
   correlation_data$pearson_pvalue[i]=as.numeric(ct_pearson$p.value)
   correlation_data$spearman_pvalue[i]=as.numeric(ct_spearman$p.value)
   correlation_data$kendall_pvalue[i]=as.numeric(ct_kendall$p.value)
   
   #data=data %>% filter(Corrected_count>0)
   #data$log_Corrected_count=log(data$Corrected_count)
   #data$log_Correct_count_scaled=data$log_Corrected_count/max(data$log_Corrected_count)
   #ggscatter(data, x = "motif_match_score", y = "log_Corrected_count", 
  #             add = "reg.line", conf.int = TRUE, 
  #             cor.coef = TRUE, cor.method = "pearson",
  #             xlab = "Motif match score", ylab = "k-mer count score")
   
   p=ggscatter(data, x = "motif_match_score_scaled", y = "Corrected_count_scaled", 
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


#3635-3453

table(!is.na(correlation_data$spearman_estimate))
#FALSE  TRUE 
#182  3453

#1. Fast & simple: plot a random subset
#You almost never need all 18M points to see the pattern.

all_data %>% filter(motif_match_score_scaled >1) %>% pull(ID) %>% unique() #None



set.seed(1)
n_sample <- 500000  # e.g. 100k points, You can tune n_sample (e.g. 50k–500k) to balance speed vs detail.
idx <- sample.int(nrow(all_data), n_sample)

plot(all_data$motif_match_score_scaled[idx], all_data$Corrected_count_scaled[idx],
     pch = ".", cex = 1,
     xlab = "x", ylab = "y")

#Use smoothScatter (base R, density-based)
#This draws a density image instead of every point.
library(viridis)
library(fields)
png(paste0(figure_path, "all_data_density_plot.png"), width = 2000, height = 2000, res = 300)
with(all_data%>% filter((motif_match_score_scaled>=-1 )& ( `Corrected_count_scaled`>=-1)), smoothScatter(motif_match_score_scaled, `Corrected_count_scaled`,
                             xlim=c(-1,1), ylim=c(-1,1),
                             nbin = 512,              # resolution of the grid, the number of equally spaced grid points for the density estimation;
                             bandwidth = 0.01,        # smoothing bandwidth, the bandwidth to be used in each coordinate direction.
                             #gridsize=128,
                             nrpoints = 0,            # how many raw points to overplot, number of points to be superimposed on the density image. The first nrpoints points from those areas of lowest regional densities will be plotted. Adding points to the plot allows for the identification of outliers. If all points are to be plotted, choose nrpoints = Inf.
                             colramp = colorRampPalette(c("white", "blue", "red")) ,
                             #colramp       = viridis, 
     #                        transformation = function(z) log1p(z),  # transform counts
                       xlab = "Motif match scaled score",
                       ylab = "k-mer count scaled score"))

dev.off()


#Hexbin / 2D histogram (aggregation)
#Using hexbin package
#install.packages("hexbin")  # if needed
library(hexbin)

# This gives a hexagon-aggregated view where color = point density.
png(paste0(figure_path, "all_data_hexbin.png"), width = 2000, height = 2000, res = 300)
with(all_data, {
  h <- hexbin(motif_match_score_scaled, `Corrected_count_scaled`, xbins = 100,
              xlab = "Motif match scaled score",
              ylab = "k-mer count scaled score" )  # adjust xbins for resolution
  plot(h)
})
dev.off()

with(all_data %>% filter((motif_match_score_scaled>=-5 )& ( `Corrected_count_scaled`>=-1)), {
  h <- hexbin(motif_match_score_scaled, `Corrected_count_scaled`, xbins = 500,
              xlab = "Motif match scaled score",
              ylab = "k-mer count scaled score" )  # adjust xbins for resolution
  plot(h)
})


ggplot(all_data%>% filter((motif_match_score_scaled>=-5 )& ( `Corrected_count_scaled`>=-1)), aes(motif_match_score_scaled, `Corrected_count_scaled`)) +
  # geom_bin2d() +  # or 
  stat_bin_2d() +
  labs( x= "Motif match scaled score",
        y = "k-mer count scaled score")

ggplot(all_data, aes(motif_match_score_scaled, `Corrected_count_scaled`)) +
  geom_binhex() +  # or stat_bin_2d()
  labs(x = "x", y = "y")

p=ggscatter(all_data%>% filter((motif_match_score_scaled>=-5 )& ( `Corrected_count_scaled`>=-1)), 
            x = "motif_match_score_scaled", y = "Corrected_count_scaled", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Scaled motif match score", ylab = "Scaled k-mer count score")


ggsave(
  filename = paste0(figure_path, "all_data.pdf"),
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)          

#Plot histogram of pearson correlations

correlation_data=correlation_data %>% filter(!is.na(pearson_estimate)) #3373 left

#correlation_data %>% filter(spearman_pvalue<0.05),
png(paste0(figure_path, "spearman_correlations.png"), width = 2000, height = 2000, res = 300)
ggplot( correlation_data, aes(x = spearman_estimate)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "", x = "Spearman correlation", y = "Frequency") +
  theme_minimal()
dev.off()


correlation_data %>% left_join(metadata %>% select(ID, study), by="ID")

ggplot( correlation_data%>% left_join(metadata %>% select(ID, study), by="ID")
, aes(x = spearman_estimate)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") + facet_wrap(~study) +
  labs(title = "", x = "Spearman correlation", y = "Frequency") +
  theme_minimal()


ggplot(correlation_data, aes(x = log10(spearman_pvalue))) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "", x = "Spearman correlation p-value", y = "Frequency") +
  theme_minimal()

min(correlation_data$spearman_pvalue) #0
max(correlation_data$spearman_pvalue) #0

#Which motifs have a high correlation

correlation_data %>% filter(spearman_estimate==max(correlation_data$spearman_estimate)) %>% pull(ID) #TGIF2_TBX21_TGGGTC40NAAG_YYI_NAGGTGTCAWN_1_3

correlation_data %>% filter(abs(spearman_estimate-0.5)<0.001) %>% pull(ID) %>% head(10) #SMAD3_HT-SELEX_TGGGTA20NGA_AC_YGTCTAGACA_1_3

correlation_data %>% filter(abs(spearman_estimate)<0.001) %>% pull(ID) %>% head() #"ETV2_FOXI1_CAP-SELEX_TACATC40NATA_AAA_RSCGGAANNRYMAACAN_1_3"

correlation_data %>% filter(spearman_estimate==min(correlation_data$spearman_estimate)) %>% pull(ID) #"ATF3_HT-SELEX_TTCATA40NTCA_KAE_NRTGACGTCAYN_1_3"


correlation_data %>% filter(abs(spearman_estimate+0.5)<0.1) %>% pull(ID) %>% head() #"ATF3_HT-SELEX_TTCATA40NTCA_KAE_NRTGACGTCAYN_1_3"
