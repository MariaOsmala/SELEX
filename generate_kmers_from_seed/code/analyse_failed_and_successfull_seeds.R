library("readr")
library("dplyr")
library("Biostrings")
library("stringdist")
library("TFBSTools")


metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

remove_ind=which(metadata$experiment=="Methyl-HT-SELEX")

failed_inds=c(999,906,920,941,950,951,992,931) #length 202

failed_inds=c(failed_inds, 
1000+c(18,81,127,136,144,145,150,153,157,174,176,177,183,193,194,205,206,209,212,
       216,220,223,224,237,238,239,240,243,244,246,247,248,250,254,255,270,271,276,
       282,285,286,287,288,289,291,292,294,295,296,297,293,281,138,
       36,46,147,152,163,165,166,170,175,178,182,185,190,192,203,207,219,229,241,256,264,290,299 
), 
1000+c(305,307,310,311,313,314,315,316,317,325,326,329,344,345,368,369,374,375,377,306,365,376, 
301,312,324,330,332,334,342,343,370,424,527),
1000+990, #Morgunava misses seed 
2000+c(221,230,254,278,279,296,76,153,286,290),
2000+c(305,306,311,376,393,424,444,557,579,308,401,502,556,559,561,589),
2000+c(607,615,616,618,635,643,646,648,650,660,665,682,684,694,695,699,705,716,718,
       719,726,733,734,749,762,772,831,832, 
      640,642,675,683,717,725,736,742),
3000+c(279,283,296,76,77,78,141),
3000+c(341,369,425,426),
3000+c(606,679,682,731,657,742,746,884),
3000+c(905,906,918)
)

successfull_inds=1:3933
successfull_inds=successfull_inds[-failed_inds]

failed_inds=failed_inds[!(failed_inds %in% remove_ind)]
successfull_inds=successfull_inds[!(successfull_inds %in% remove_ind)]


failed_seeds=metadata$seed[failed_inds] #197

successfull_seeds=metadata$seed[successfull_inds] #3439

length(failed_seeds) #197
length(successfull_seeds) #3439
#sum 3636

write.table(metadata$ID[successfull_inds], 
            "/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/successfull_exps_excluding_Methyl-HT-SELEx.txt",
            quote=FALSE, 
            sep="\t",
            row.names=FALSE,
            col.names=FALSE)


hist(sapply(failed_seeds, nchar))

hist(sapply(successfull_seeds, nchar))


min(nchar(failed_seeds), na.rm=TRUE ) #18
max(nchar(failed_seeds), na.rm=TRUE ) #32

min(nchar(successfull_seeds), na.rm=TRUE ) #7
max(nchar(successfull_seeds), na.rm=TRUE ) #31

count_degenerate <- function(seq) {
  # seq: character string (k-mer)
  nchar(gsub("[ACGT]", "", seq))
}

count_N <- function(seq) {
  nchar(gsub("[^Nn]", "", seq))   # keep only N/n, count length
}

# Example

hist(sapply(failed_seeds, count_degenerate))
hist(sapply(successfull_seeds, count_degenerate))

min(sapply(failed_seeds, count_degenerate), na.rm=TRUE) #9 
max(sapply(failed_seeds, count_degenerate), na.rm=TRUE) #22

min(sapply(successfull_seeds, count_degenerate), na.rm=TRUE) #0
max(sapply(successfull_seeds, count_degenerate), na.rm=TRUE) #21


hist(sapply(failed_seeds, count_N))
hist(sapply(successfull_seeds, count_N))

min(sapply(failed_seeds, count_N), na.rm=TRUE) #5
max(sapply(failed_seeds, count_N), na.rm=TRUE) #19

min(sapply(successfull_seeds, count_N), na.rm=TRUE) #0
max(sapply(successfull_seeds, count_N), na.rm=TRUE) #19

