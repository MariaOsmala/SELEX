library(readr)
library(dplyr)
library(stringr)
data <- read_delim("/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv", 
                   delim = "\t", 
                   escape_double = FALSE, 
                   trim_ws = TRUE)



data=data %>% filter(study=="Jolma2015")

included_signal=do.call(rbind, str_split(data$CSC_SELEX_filename, "/"))[,8] %>% unique() #327

included_background=as.data.frame(do.call(rbind, str_split(data$CSC_SELEX_background_filename, "/")))

included_background=included_background %>% filter(V6=="Jolma2015")
included_background=included_background[,8] %>% unique

included=c(included_signal, included_background)

write.table(included, file="/projappl/project_2013895/SELEX/data-to-allas/Jolma2015.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


data=data %>% filter(study=="Yin2017")

included_signal=do.call(rbind, str_split(data$CSC_SELEX_filename, "/"))[,8] %>% unique() 

included_background=as.data.frame(do.call(rbind, str_split(data$CSC_SELEX_background_filename, "/")))

included_background=included_background %>% filter(V6=="Yin2017")
included_background=included_background[,8] %>% unique

included=c(included_signal, included_background)

write.table(included, file="/projappl/project_2013895/SELEX/data-to-allas/Yin2017.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


#all

urls=unique(c(data$Allas_SELEX_filename, data$Allas_SELEX_background_filename))

write.table(urls, file="/projappl/project_2013895/SELEX/data-to-allas/all.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
