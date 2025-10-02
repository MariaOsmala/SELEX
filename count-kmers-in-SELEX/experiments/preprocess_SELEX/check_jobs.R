library("readr")
library("dplyr")
library("tidyverse")


metadata=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv", delim="\t")

check_jobs=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/check_jobs.tsv", delim="\t")


check_jobs=check_jobs%>%
  separate(JobID, into = c("array_job", "array_index"), sep = "_", remove = FALSE)

check_jobs <- check_jobs %>%
  mutate(across(c( "array_index"  ),
                as.numeric))


# check_jobs <- check_jobs %>%
#   mutate(across( "array_job",~factor(.x, levels = c("29343847","29345862","29345938","29347707","29359880",
#                                                     "29347902","29347983","29348008","29348079","29355999",
#                                                     "29356094","29356233","29357154","29357891","29357964","29358009"))
#   ))
# 

check_jobs=check_jobs %>%
  arrange(array_job, array_index)

check_jobs <- check_jobs %>%
  mutate(rownum = row_number()) %>%
  relocate(rownum, .before = 1)  

check_jobs$State %>% table()
#CANCELLED COMPLETED    FAILED 
#4       274        80 

steps0=read_delim("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/steps0.csv", delim="|")

steps0=steps0%>%
  separate(JobID, into = c("array_job", "array_index"), sep = "_", remove = FALSE)

steps0 <- steps0 %>%
  mutate(across(c( "array_index"  ),
                as.numeric))


# steps0 <- steps0 %>%
#   mutate(across( "array_job",~factor(.x, levels = c("29343847","29345862","29345938","29347707","29359880",
#                                                     "29347902","29347983","29348008","29348079","29355999",
#                                                     "29356094","29356233","29357154","29357891","29357964","29358009"))
#   ))
# 

steps0=steps0 %>%
  arrange(array_job, array_index)

steps0 <- steps0 %>%
  mutate(rownum = row_number()) %>%
  relocate(rownum, .before = 1)  


check_jobs$State %>% table()
#CANCELLED COMPLETED    FAILED 
#4       274        80 

steps0$State %>% table()
#CANCELLED CANCELLED by 0      COMPLETED         FAILED 
#4             56            218             80 

check_jobs=check_jobs %>% mutate(State_step0=steps0$State,.after = State)
