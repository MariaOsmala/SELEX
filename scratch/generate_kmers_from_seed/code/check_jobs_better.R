library("dplyr")
library("tidyverse")
rm(list=ls())

df=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/jobs_summary.tsv")
df=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/failed.tsv")

df=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/kmers_and_Hamming12.tsv")
df=df %>%filter(!(State=="PENDING"))

df=read_delim("/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/Hamming12.tsv")

df=df%>%
  separate(JobID, into = c("array_job", "array_index"), sep = "_", remove = FALSE)

df <- df %>%
  mutate(across(c( "array_index"  ),
                as.numeric))


df <- df %>%
  mutate(across( "array_job",~factor(.x, levels = c("29343847","29345862","29345938","29347707","29359880",
                                                    "29347902","29347983","29348008","29348079","29355999",
                                                    "29356094","29356233","29357154","29357891","29357964","29358009"))
                ))




df=df %>%
  arrange(array_job, array_index)

df <- df %>%
  mutate(rownum = row_number()) %>%
  relocate(rownum, .before = 1)  


metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

metadata=metadata %>% filter(!(experiment %in% "Methyl-HT-SELEX")) #3636

df$motif=metadata$ID[1:nrow(df)]
df$seed=metadata$seed[1:nrow(df)]

df <- df %>% relocate(motif, .after=rownum)
df <- df %>% relocate(seed, .after=motif)

df <- df %>%
  rename("ReqMem_K"="ReqMem")

df <- df %>%
  rename("MaxRSS_K"="MaxRSS")

df$ReqMem_K=gsub("K", "", df$ReqMem_K)
df$MaxRSS_K=gsub("K", "", df$MaxRSS_K)

df <- df %>%
  mutate(across(c( "ReqMem_K", "MaxRSS_K"  ),
                as.numeric))


df$State %>% table()
#CANCELLED COMPLETED    FAILED   RUNNING 
#89      3455        31        61 

# Hamming 1 k-mers 
# CANCELLED means time runts out
# CANCELLED COMPLETED    FAILED   
# Mem 10G is enough, could be increased to 15
#CANCELLED COMPLETED    FAILED 
#149      3456        31 


df$State[which(df$rownum%in% c(21,63,180,294,297,300,302,306,384,386,409,411,619,658,685,686,706,730,738,819,830,832,841,842,869,892,893))]="COMPLETED"

View(df %>% filter(State%in% c("FAILED", "CANCELLED"))) 

df %>% filter(State%in% c( "CANCELLED")) %>% pull(seed)
#took max 2-22:40:04

write.table(df, 
            "/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/failed_experiments_Hamming1.tsv",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)



df %>% filter(State%in% c("FAILED", "CANCELLED")) %>% pull(motif)

# COMPLETED    FAILED   RUNNING 
# 3344        28       264 


df %>% filter(State=="FAILED") %>% pull(rownum)



paste(df %>% filter(State == "FAILED") %>% pull(rownum), collapse = ",")
paste(df %>% filter(State == "FAILED") %>% pull(array_index), collapse = ",")
