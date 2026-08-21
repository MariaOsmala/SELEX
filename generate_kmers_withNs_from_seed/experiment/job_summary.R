library("readr")
library("dplyr")
library("tidyverse")


check_jobs=read_delim("/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/experiment/jobs_summary.tsv", delim="\t")


check_jobs=check_jobs%>%
  separate(JobID, into = c("array_job", "array_index"), sep = "_", remove = FALSE)

check_jobs <- check_jobs %>%
  mutate(across(c( "array_index"  ),
                as.numeric))

check_jobs=check_jobs %>%
  arrange(array_job, array_index)

check_jobs <- check_jobs %>%
  mutate(rownum = row_number()) %>%
  relocate(rownum, .before = 1)  

check_jobs$State %>% table()
#CANCELLED COMPLETED    FAILED 
#4       274        80 

steps0=read_delim("/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/experiment/jobs_summary_steps0.csv", delim="|")

steps0=steps0%>%
  separate(JobID, into = c("array_job", "array_index"), sep = "_", remove = FALSE)

steps0 <- steps0 %>%
  mutate(across(c( "array_index"  ),
                as.numeric))

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

check_jobs$State_step0 %>% table()

check_jobs %>% filter(State_step0=="CANCELLED by 0") %>% nrow() #106

paste(check_jobs %>% filter(State_step0=="CANCELLED by 0") %>% pull(array_index), collapse=",")


collapse_ranges <- function(x, sep = ",", dash = "-", trailing_comma = TRUE) {
  v <- sort(unique(na.omit(as.integer(x))))
  if (length(v) == 0) return(if (trailing_comma) "" else character())
  # find run starts/ends
  starts <- c(1, which(diff(v) != 1) + 1)
  ends   <- c(starts[-1] - 1, length(v))
  parts <- mapply(function(a, b) if (a == b) as.character(v[a]) else paste0(v[a], dash, v[b]),
                  starts, ends, USE.NAMES = FALSE)
  out <- paste(parts, collapse = paste0(sep))
  if (trailing_comma) paste0(out, sep) else out
}

collapse_ranges(check_jobs %>% filter(State_step0=="CANCELLED by 0") %>% pull(array_index))


238,395-396,398,412,454,467,469,                    472-478,482,486,489-490,492,    512-514,517-518,522-523,526,531-532,535-536,                    539-545,            547-551,        554-557,559,    561-563,                                                                   567-585,         587-590,    592-594,                597-602,604,606-607,                            609-635,"

238,395,396,398,412,454,467,469,472,473,474,475,476,477,478,482,486,489,490,492,512,513,514,517,518,522,523,526,531,532,535,536,539,540,541,542,543,544,545,547,548,549,550,551,554,555,556,557,559,561,562,563,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,587,588,589,590,592,593,594,597,598,599,600,601,602,604,606,607,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635"

