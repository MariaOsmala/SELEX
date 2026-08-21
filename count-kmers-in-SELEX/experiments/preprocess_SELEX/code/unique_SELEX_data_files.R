library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# metadata_all=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026.tsv", delim="\t") #3933
# metadata_all=metadata_all %>% filter((experiment %in% c("HT-SELEX", "CAP-SELEX")) &
#                           study!="Morgunova2015")  #3635
# 
# metadata_all %>% filter(CSC_SELEX_filename!="" & CSC_SELEX_background_filename!="") %>% nrow() #3596

metadata=read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_19032026.tsv", delim="\t") #3596

all_files=unique(c(metadata$CSC_SELEX_filename, metadata$CSC_SELEX_background_filename)) #4552

#Are the filenames the same
file_names=gsub(".fastq.gz", "", do.call(rbind, strsplit(all_files, split="/"))[,8])

unique(file_names) %>% length() #4552

write.table(all_files, file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/unique_SELEX.txt", 
            col.names = FALSE, row.names = FALSE, quote=FALSE)


# Convert strings like "10 Mbp" or "850 Kbp" to numeric bases
parse_total_bases <- function(x) {
  x <- trimws(x)
  num  <- as.numeric(sub("^([0-9.]+).*$", "\\1", x))
  unit <- sub("^[0-9.]+\\s*", "", x)
  
  mult <- c("bp"  = 1,
    "kbp" = 1e3,
    "Mbp" = 1e6,
    "Gbp" = 1e9
  )
  
  num * unname(mult[unit])
}

read_fastqc_basic <- function(fastqc_dir) {
  #fastqc_dir=qc_dirs[1]
  f <- file.path(fastqc_dir, "fastqc_data.txt")
  x <- readLines(f, warn = FALSE)
  
  start <- grep("^>>Basic Statistics", x)
  end   <- grep("^>>END_MODULE", x)
  end   <- end[end > start][1]
  
  block <- x[(start + 1):(end - 1)]
  
  tab <- read.delim(
    text = paste(block, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(tab)=c("Measure", "Value")
  vals <- setNames(tab$Value, tab$Measure)
  
  data.frame(
    sample = sub("_fastqc$", "", basename(fastqc_dir)),
    total_sequences = as.numeric(vals[["Total Sequences"]]),
    total_bases_raw = vals[["Total Bases"]],
    total_bases = parse_total_bases(vals[["Total Bases"]]),
    poor_quality_sequences = as.numeric(vals[["Sequences flagged as poor quality"]]),
    sequence_length = vals[["Sequence length"]],
    pct_gc = as.numeric(vals[["%GC"]]),
    stringsAsFactors = FALSE
  )
}

read_fastqc_summary <- function(fastqc_dir) {
  f <- file.path(fastqc_dir, "summary.txt")
  
  df=read.delim(
    f,
    sep = "\t",
    header = FALSE,
    col.names = c("status", "module", "filename"),
    stringsAsFactors = FALSE
  ) |>
    mutate(sample = sub("_fastqc$", "", basename(fastqc_dir))) |>
    select(sample, status, module, filename)
  
  df %>%
    select(sample, module, status) %>%
    pivot_wider(
      names_from = module,
      values_from = status
    )
  
}

# Top folder containing many *_fastqc directories
qc_root <- "/scratch/project_2013895/SELEX/fastqc"

qc_dirs <- list.dirs(qc_root, recursive = FALSE, full.names = TRUE) #4543
qc_dirs <- qc_dirs[grepl("_fastqc$", qc_dirs)]

#data with results
successfull=gsub("_fastqc", "", do.call(rbind, strsplit(qc_dirs, split="/"))[,6]) #4543, 9 fastqs are corrupted
file_names[which(!(file_names %in% successfull))]

#The fastq file is corrupted
# "ELF1_CREM_3_YLIII_TTCATT40NTAT"  1747
# "ELF4_CREM_3_YLIII_TCTAGA40NCGA"  1768
# "HOXA4_CREM_3_YLIII_TAGACT40NGTG" 1987
# "HOXD3_CREM_3_YLIII_TGTAAG40NACT" 2074 
# "HOXD4_CREM_3_YLIII_TTCTGG40NTGC" 2080 
# "PAX1_VSX2_3_YLIIII_TATAGC40NGAC" 2263 
# "PAX2_VSX2_3_YLIIII_TTGCTG40NCTA" 2283 
# "POU4F3_CREM_3_YLIII_TAAGTT40NGAT" 2346
# "RFX5_CREM_3_YLIII_TCCCTT40NCAT"   2367




basic_info <- bind_rows(lapply(qc_dirs, read_fastqc_basic))
summary_info <- bind_rows(lapply(qc_dirs, read_fastqc_summary))

df=data.frame(data_ID=file_names, filenames=all_files)

df$type=NA

df$type[(df$filenames %in% metadata$CSC_SELEX_filename) & !(df$filenames %in% metadata$CSC_SELEX_background_filename)]="signal"
df$type[!(df$filenames %in% metadata$CSC_SELEX_filename) & (df$filenames %in% metadata$CSC_SELEX_background_filename)]="background"
df$type[(df$filenames %in% metadata$CSC_SELEX_filename) & (df$filenames %in% metadata$CSC_SELEX_background_filename)]="both"

df$type %>% table(useNA="always")
#background       both     signal       <NA> 
#  1878         61       2613          0 


df=df %>% left_join(basic_info, by=c("data_ID"="sample"))
df=df %>% left_join(summary_info, by=c("data_ID"="sample"))


write_delim(
  df,
  file="/projappl/project_2013895/SELEX/quality_control/fastqc/quality_control.tsv",
  delim = "\t"
)


summary(df)

names(df)

is.na(df$total_sequences) %>% table()
#FALSE  TRUE 
#4543     9 (corrupted fastq files, these are all signal files)

df$sequence_length=as.numeric(df$sequence_length)


plot_histograms <- function(data, cols, bins = 30, scales = "fixed") {
  data_long <- data %>%
    select(all_of(cols)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  ggplot(data_long, aes(x = value)) +
    geom_histogram(bins = bins, color = "black", fill = "grey80") +
    facet_wrap(~ variable, scales = "free") +
    scale_x_continuous(labels = label_comma()) +
    theme_bw()
}


#are these any NA values

is.na(df$total_sequences) %>% table() #9
is.na(df$total_bases) %>% table() #9
is.na(df$sequence_length) %>% table() #9
is.na(df$pct_gc) %>% table() #9

p=plot_histograms(df, cols=c("total_sequences","total_bases", 
                           #"poor_quality_sequences", #0 poor quality sequences
                           "sequence_length", "pct_gc"))

ggsave("/projappl/project_2013895/SELEX/quality_control/fastqc/fastqc_numerical.pdf", plot = p, width = 8, height = 6)


plot_categorical_bars <- function(data, cols) {
  data_long <- data %>%
    select(all_of(cols)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "value"
    )  %>%
    mutate(value = ifelse(is.na(value), "Missing", as.character(value)))
  
  ggplot(data_long, aes(x = value)) +
    geom_bar() +
    facet_wrap(~ variable, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

df$`Basic Statistics` %>% table(useNA="always") #9NAs
df$`Per base sequence quality` %>% table(useNA="always") #9NAs
df$`Per tile sequence quality` %>% table(useNA="always") #1114 NAs
df$`Per sequence quality scores` %>% table(useNA="always") #9 NAs
df$`Per base sequence content` %>% table(useNA="always") #9 NAs
df$`Per sequence GC content`%>% table(useNA="always") #9 NAs
df$`Per base N content` %>% table(useNA="always") #9 NAs
df$`Sequence Length Distribution` %>% table(useNA="always") #9 NAs
df$`Sequence Duplication Levels` %>% table(useNA="always") #9 NAs
df$`Overrepresented sequences`%>% table(useNA="always") #9 NAs
df$`Adapter Content` %>% table(useNA="always") #9 NAs

p=plot_categorical_bars(df, cols=c(
                                 #"Basic Statistics", #All pass
                                 #"Per base sequence quality", All pass
                                 "Per tile sequence quality",   
                                 #"Per sequence quality scores", All pass
                                 "Per base sequence content",
                                 "Per sequence GC content",   
                                 "Per base N content",          
                                 #"Sequence Length Distribution", #All pass
                                 "Sequence Duplication Levels",
                                 "Overrepresented sequences"
                                 #"Adapter Content" #All pass
                                 )
                                 )  

ggsave("/projappl/project_2013895/SELEX/quality_control/fastqc/fastqc_categorical.pdf", plot = p, width = 8, height = 6)


df %>% filter(`Per base N content`=="WARN") %>% pull(data_ID)


# Quality measures after removal of reads with N --------------------------


qc_root <- "/scratch/project_2013895/SELEX/fastqc_reads_with_N_removed"

qc_dirs <- list.dirs(qc_root, recursive = FALSE, full.names = TRUE) #4544
qc_dirs <- qc_dirs[grepl("_fastqc$", qc_dirs)]

#data with results
successfull=gsub("_fastqc", "", do.call(rbind, strsplit(qc_dirs, split="/"))[,6]) #4543, 8 fastqs are corrupted, one less than above
file_names[which(!(file_names %in% successfull))]


basic_info <- bind_rows(lapply(qc_dirs, read_fastqc_basic))
summary_info <- bind_rows(lapply(qc_dirs, read_fastqc_summary))

df_after_filtering=data.frame(data_ID=file_names, filenames=all_files)



df_after_filtering= df_after_filtering%>% left_join(basic_info, by=c("data_ID"="sample"))
df_after_filtering=df_after_filtering %>% left_join(summary_info, by=c("data_ID"="sample"))

names(df_after_filtering)[-c(1,2)]=paste0(names(df_after_filtering)[-c(1,2)], " after filtering")

df_after_filtering$filenames=NULL

df = df %>% left_join(df_after_filtering, by="data_ID")

write_delim(
  df,
  file="/projappl/project_2013895/SELEX/quality_control/fastqc/quality_control_after_filtering.tsv",
  delim = "\t"
)




names(df)

is.na(df$`total_sequences after filtering`) %>% table()
#FALSE  TRUE 
#4544     8 (corrupted fastq files, these are all signal files)

df$`sequence_length after filtering`=as.numeric(df$`sequence_length after filtering`)




p=plot_histograms(df, cols=c("total_sequences after filtering","total_bases after filtering", 
                            # "poor_quality_sequences after filtering", #0 poor quality sequences
                             "sequence_length after filtering", "pct_gc after filtering"))

ggsave("/projappl/project_2013895/SELEX/quality_control/fastqc/fastqc_numerical_after_filtering.pdf", plot = p, width = 8, height = 6)




p=plot_categorical_bars(df, cols=c(
  #"Basic Statistics after filtering", #All pass
  #"Per base sequence quality after filtering", #All pass
  "Per tile sequence quality after filtering",   
 # "Per sequence quality scores after filtering", #All pass
  "Per base sequence content after filtering",
  "Per sequence GC content after filtering",   
  "Per base N content after filtering",          
 # "Sequence Length Distribution after filtering", #All pass
  "Sequence Duplication Levels after filtering",
  "Overrepresented sequences after filtering"
  #"Adapter Content after filtering" #All pass
)
)  

ggsave("/projappl/project_2013895/SELEX/quality_control/fastqc/fastqc_categorical_after_filtering.pdf", plot = p, width = 8, height = 6)


df %>% filter(`Overrepresented sequences after filtering`=="FAIL") %>% pull(data_ID)

#For this data unique reads are used

uniq_exps=metadata %>% filter(unique_background==TRUE) %>% pull(CSC_SELEX_filename, CSC_SELEX_background_filename) %>% unique()

test=df %>% filter(filenames %in% uniq_exps) %>% select(`Overrepresented sequences`, 
                                                        `Overrepresented sequences after filtering`, 
                                                        `Sequence Duplication Levels`, 
                                                        `Sequence Duplication Levels after filtering`, 
                                                        data_ID)

#For which motifs the lambda computation failed, is it because of issues in the data?



lambda_info=read_delim("/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_lambda_info.tsv", delim="\t")

metadata=metadata %>% left_join(lambda_info %>% select(ID, lambda, N_sig, N_bg), by="ID")

metadata=metadata %>% left_join(df %>% select(filenames,   "total_sequences after filtering",             
                                              "total_bases_raw after filtering",
                                              "total_bases after filtering",
                                              "poor_quality_sequences after filtering",      
                                              "sequence_length after filtering",
                                              "pct_gc after filtering",
                                              "Basic Statistics after filtering",
                                              "Per base sequence quality after filtering",
                                              "Per tile sequence quality after filtering",
                                              "Per sequence quality scores after filtering", 
                                              "Per base sequence content after filtering",
                                              "Per sequence GC content after filtering",
                                              "Per base N content after filtering",          
                                              "Sequence Length Distribution after filtering",
                                              "Sequence Duplication Levels after filtering",
                                              "Overrepresented sequences after filtering",   
                                              "Adapter Content after filtering") %>%rename_with(~ paste0(.x, "_signal")),
                                by=c("CSC_SELEX_filename"="filenames_signal")
)


metadata=metadata %>% left_join(df %>% select(filenames,   "total_sequences after filtering",             
                                              "total_bases_raw after filtering",
                                              "total_bases after filtering",
                                              "poor_quality_sequences after filtering",      
                                              "sequence_length after filtering",
                                              "pct_gc after filtering",
                                              "Basic Statistics after filtering",
                                              "Per base sequence quality after filtering",
                                              "Per tile sequence quality after filtering",
                                              "Per sequence quality scores after filtering", 
                                              "Per base sequence content after filtering",
                                              "Per sequence GC content after filtering",
                                              "Per base N content after filtering",          
                                              "Sequence Length Distribution after filtering",
                                              "Sequence Duplication Levels after filtering",
                                              "Overrepresented sequences after filtering",   
                                              "Adapter Content after filtering") %>%
                                  rename_with(~ paste0(.x, "_background")),
                                by=c("CSC_SELEX_background_filename"="filenames_background")
)


test=metadata %>% filter(is.na(lambda)) %>% select(lambda, N_sig, N_bg, 

