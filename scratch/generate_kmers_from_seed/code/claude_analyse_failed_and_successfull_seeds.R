library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")
library("ggplot2")

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t") #this contains all motifs plus SELEX data info
metadata= metadata %>% filter(experiment %in% c("HT-SELEX", "CAP-SELEX"))
metadata=metadata %>% filter(!is.na(seed)) #3635

#sort based on the seed length
metadata=metadata[order(metadata$length),]


results_path="/scratch/project_2013895/SELEX/combined_kmers_from_long_degenerate_seeds"

files <- list.files(results_path, pattern = ".txt", full.names = TRUE) #3544


get_n <- function(f) {
  out <- system(sprintf("wc -l -- %s", shQuote(f)), intern = TRUE)
  as.integer(sub("^\\s*(\\d+).*", "\\1", out))  # take the count only
}

# jellyfish succeeded for 6966/2=3483

n <- vapply(files, get_n, integer(1))

names(n)=gsub("/scratch/project_2013895/SELEX/combined_kmers_from_long_degenerate_seeds/", "", names(n))
names(n)=gsub(".txt", "", names(n))
row_numbers=data.frame(ID=names(n), number=as.integer(n))


metadata$kmer_generation_successfull=TRUE

metadata$kmer_generation_successfull[metadata$ID %in% (metadata %>% filter(!(ID %in% row_numbers$ID)) %>% pull(ID))]=FALSE

metadata %>% select(kmer_generation_successfull) %>% table()
#FALSE  TRUE 
#91  3544 

View(metadata %>% select(ID, seed, length,kmer_generation_successfull))

#For the successfull k-mer generation, how many k-mers were generated

metadata$kmer_nro=NA
match_ind=match(row_numbers$ID,metadata$ID)

#head(row_numbers$ID)
#head(metadata$ID[match_ind])

metadata$kmer_nro[match_ind]=row_numbers$number

write.table(metadata, 
            "/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/claude_Hamming1_kmercount_analysis.tsv",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)




info=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/Hamming1_score_analysis.tsv")
claude_info=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/claude_Hamming1_kmercount_analysis.tsv")

info %>% nrow() #3461


claude_info$kmer_generation_successfull %>% table()

#FALSE  TRUE 
#91  3544

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

degeneracy_stats_list=lapply(claude_info$seed, degeneracy_stats)               




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

for(i in 1:nrow(claude_info)){
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




claude_info=left_join(claude_info, degeneracy_stats_df, by="ID")

p=ggplot(claude_info %>% filter(kmer_generation_successfull==TRUE), aes(x = length, y = kmer_nro, 
                                                                   color=D_norm)) +
  scale_color_gradient(low = "blue", high = "red")+
  geom_jitter(width = 0.2, height = 0.2)+
  scale_y_log10() +
  labs(title = "Hamming 1 distance k-mers", x = "seed length", y = "k-mer number") +
  theme_minimal()




p=ggplot(claude_info %>% filter(kmer_generation_successfull==FALSE), aes(x = length, y = D_norm) )+
  #scale_color_gradient(low = "blue", high = "red")+
  geom_point()+xlim(c(18,26))+ylim(c(0.35, 0.65))+
  #geom_jitter(width = 0.2, height = 0.2)+
  #scale_y_log10() +
  labs(title = "Hamming 1 distance k-mers", x = "seed length", y = "D_norm") +
  theme_minimal()


ggsave(
  filename = "/projappl/project_2013895/SELEX/generate_kmers_from_seed/Figures/claude_Hamming1_failed_seeds_seedlength_vs_normalized_degeneracy.pdf",
  plot = p,
  device = cairo_pdf,      # nicer text/antialiasing; needs Cairo installed
  width = 7, height = 5, units = "in"  # set page size
)    

#Shortest k-mers for which the k-mer generation failed

write.table(claude_info, 
            "/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/claude_metadata_failed_and_successfull_seeds_degeneracy.tsv",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)


claude_info=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/claude_metadata_failed_and_successfull_seeds_degeneracy.tsv")

claude_info %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(seed)
# "NNNACGANNNNNNTCGTNNN"

# Earlier 
#"NNCCGGNNNNNNCCGGNN"
