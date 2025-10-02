library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("ggplot2")

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

metadata=metadata %>% filter(!(experiment %in% "Methyl-HT-SELEX")) 

is.na(metadata$seed) %>% table()

metadata=metadata %>% filter(!(ID=="E2F8_Morgunova2015")) #3635

results=dir("/scratch/project_2013895/SELEX/streamed_kmers_exclude_Methyl-HT-SELEX/")
results=gsub(".tsv", "", results) #3461

metadata$kmer_generation_successfull=TRUE

metadata$kmer_generation_successfull[metadata$ID %in% (metadata %>% filter(!(ID %in% results)) %>% pull(ID))]=FALSE

metadata %>% select(kmer_generation_successfull) %>% table()
#FALSE  TRUE 
#174  3461 

#For the successfull k-mer generation, how many k-mers were generated

successfull=metadata %>% filter(kmer_generation_successfull==TRUE)

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



info=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming1_score_analysis.tsv")







ggplot(info, aes(x = kmer_nro)) +
  geom_histogram(bins=100) + 
  scale_x_continuous() + #xlim(c(0,1000))+ 
  scale_x_log10()+ #scale_x_log10()+
  labs(
    title = "Hamming1",
    x = "k-mer number",
    y = "count"
  ) +
  theme_minimal() +theme(plot.title = element_text(size = 10) )






info=info %>% left_join(metadata %>% select(ID, length, IC))
               
         
p=ggplot(info, aes(x = length, y = kmer_nro, color=IC/length)) +
      scale_color_gradient(low = "blue", high = "red")+
      geom_jitter(width = 0.2, height = 0.2)+
      scale_y_log10() +
      labs(title = "Hamming 1 distance k-mers", x = "seed length", y = "k-mer number") +
      theme_minimal()
                        
ggsave(
  filename = "/projappl/project_2013895/SELEX/generate_kmers_from_seed/Figures/Hamming1_kmer_number_vs_seedlength_averagemotifIC.pdf",
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)                       
               
#Compute degeneracy scores for seeds:                         
                        
# IUPAC -> set of concrete bases
.iupac_map <- list(
  A=c("A"), C=c("C"), G=c("G"), T=c("T"),
  R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
  K=c("G","T"), M=c("A","C"),
  B=c("C","G","T"), D=c("A","G","T"), H=c("A","C","T"), V=c("A","C","G"),
  N=c("A","C","G","T")
)
.comp_base <- c(A="T", C="G", G="C", T="A")

parse_iupac_sets <- function(pat) {
  chars <- strsplit(toupper(pat), "", fixed = TRUE)[[1]]
  sets <- lapply(seq_along(chars), function(i) {
    s <- .iupac_map[[chars[i]]]
    if (is.null(s)) stop(sprintf("Unsupported IUPAC symbol '%s' at pos %d", chars[i], i))
    s
  })
  sets
}

degeneracy_stats <- function(pat, big = TRUE) {
  sets <- parse_iupac_sets(pat)
  k <- length(sets)
  sizes <- vapply(sets, length, integer(1))
  
  # Cardinality (optionally big integer)
  if (big) {
    stopifnot(requireNamespace("gmp", quietly = TRUE))
    C <- Reduce(`*`, lapply(sizes, gmp::as.bigz), init = gmp::as.bigz(1))
    log2C <- as.numeric(sum(log2(sizes)))  # may be Inf for very large k
  } else {
    C <- prod(sizes)
    log2C <- sum(log2(sizes))
  }
  
  # Normalized degeneracy in [0,1]
  D_norm <- mean(log(sizes, base = 4))
  
  # Gini-style impurity (normalized)
  G_norm <- mean((1 - 1/sizes) / (1 - 1/4))
  
  # RC-canonical count via Burnside
  if (k %% 2L == 1L) {
    F <- if (big) gmp::as.bigz(0) else 0L
  } else {
    pair_counts <- vapply(seq_len(k/2), function(i) {
      Si <- sets[[i]]
      Sj <- sets[[k + 1L - i]]
      compSj <- unname(.comp_base[Sj])
      length(intersect(Si, compSj))
    }, integer(1))
    if (big) {
      F <- Reduce(`*`, lapply(pair_counts, gmp::as.bigz), init = gmp::as.bigz(1))
    } else {
      F <- prod(pair_counts)
    }
  }
  
  C_canon <- if (k %% 2L) {
    if (big) C / gmp::as.bigz(2) else C / 2
  } else {
    if (big) (C + F) / gmp::as.bigz(2) else (C + F) / 2
  }
  
  list(
    k = k,
    per_position_sizes = sizes,
    count = C,
    bits = log2C,
    D_norm = D_norm,
    G_norm = G_norm,
    canonical_count = C_canon
  )
}


degeneracy_stats("ATGC", big=FALSE)         # fully specific
degeneracy_stats("NNNN")         # fully degenerate
degeneracy_stats("ATGCNRY")      # mixed
degeneracy_stats("ACGT", big=TRUE) # big-int safe

degeneracy_stats_list=lapply(metadata$seed, degeneracy_stats)               

degeneracy_stats_df=tibble(ID=rep(NA, length(degeneracy_stats_list)), 
                               k=rep(NA, length(degeneracy_stats_list)),
                               per_position_sizes=rep("", length(degeneracy_stats_list)),
                               count=rep(gmp::as.bigz(NA), length(degeneracy_stats_list)),
                               bits=rep(NA, length(degeneracy_stats_list)),
                               D_norm=rep(NA, length(degeneracy_stats_list)),
                               G_norm=rep(NA, length(degeneracy_stats_list)) #,
                               canonical_count=rep(gmp::as.bigz(NA), length(degeneracy_stats_list))
                           ) 


n <- length(degeneracy_stats_list)

degeneracy_stats_df <- tibble::tibble(
  ID = NA_character_,
  k  = NA_integer_,
  # this one should be a list-column if you plan to store integer vectors:
  per_position_sizes = NA_character_,
  count            = NA_character_,
  bits             = NA_real_,
  D_norm           = NA_real_,
  G_norm           = NA_real_,
  canonical_count = NA_character_,
)

degeneracy_stats_df=degeneracy_stats_df[1,]

degeneracy_stats_df=degeneracy_stats_df[rep(1, n), ]

for(i in 1:nrow(metadata)){
  print(i)
  degeneracy_stats_df$ID[i]=metadata$ID[i]
  degeneracy_stats_df$k[i]=degeneracy_stats_list[[i]]$k
  degeneracy_stats_df$per_position_sizes[i]=paste0(as.character(degeneracy_stats_list[[1]]$per_position_sizes), collapse=",")
  degeneracy_stats_df$count[i]=as.character(degeneracy_stats_list[[i]]$count)
  degeneracy_stats_df$bits[i]=degeneracy_stats_list[[i]]$bits
  degeneracy_stats_df$D_norm[i]=degeneracy_stats_list[[i]]$D_norm
  degeneracy_stats_df$G_norm[i]=degeneracy_stats_list[[i]]$G_norm
  degeneracy_stats_df$canonical_count[i]=as.character(degeneracy_stats_list[[i]]$canonical_count)
  
}  


metadata= left_join(metadata, info %>% select(ID,kmer_nro, min_score, max_score, min_scaled_score, max_scaled_score), by="ID")

metadata=left_join(metadata, degeneracy_stats_df, by="ID")

p=ggplot(metadata %>% filter(kmer_generation_successfull==TRUE), aes(x = length, y = kmer_nro, 
                                                                   color=D_norm)) +
  scale_color_gradient(low = "blue", high = "red")+
  geom_jitter(width = 0.2, height = 0.2)+
  scale_y_log10() +
  labs(title = "Hamming 1 distance k-mers", x = "seed length", y = "k-mer number") +
  theme_minimal()


ggsave(
  filename = "/projappl/project_2013895/SELEX/generate_kmers_from_seed/Figures/Hamming1_kmer_number_vs_seedlength_normalized_degeneracy.pdf",
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)           


p=ggplot(metadata %>% filter(kmer_generation_successfull==FALSE), aes(x = length, y = D_norm) )+
  #scale_color_gradient(low = "blue", high = "red")+
  geom_point()+
  #geom_jitter(width = 0.2, height = 0.2)+
  #scale_y_log10() +
  labs(title = "Hamming 1 distance k-mers", x = "seed length", y = "D_norm") +
  theme_minimal()


ggsave(
  filename = "/projappl/project_2013895/SELEX/generate_kmers_from_seed/Figures/Hamming1_failed_seeds_seedlength_vs_normalized_degeneracy.pdf",
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)    

