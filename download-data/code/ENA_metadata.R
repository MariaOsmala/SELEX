library("readr")
library("dplyr")
metadata=read_delim("/scratch/project_2013895/SELEX/data/Yin2017/filereport_read_run_PRJEB9797_tsv.txt")


metadata[1:5,c("run_accession","library_strategy","library_source","library_selection","scientific_name","instrument_model","sample_alias",  "sample_title" )]


metadata$library_strategy %>% table()
metadata$library_source %>% table()
metadata$library_selection %>% table()
metadata$scientific_name %>% table()
metadata$instrument_model %>% table()
