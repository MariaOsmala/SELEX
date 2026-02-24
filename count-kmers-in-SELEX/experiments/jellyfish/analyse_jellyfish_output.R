library("readr")
library("dplyr")
#4^7=16384
#4^7/2=8192

signal=read_delim("/scratch/project_2013895/SELEX/jellyfish/MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6_signal_mer_counts_dumps.txt", col_names = FALSE) #8192
background=read_delim("/scratch/project_2013895/SELEX/jellyfish/MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6_background_mer_counts_dumps.txt", col_names=FALSE) #8192

kmer_signal=read_delim("/scratch/project_2013895/SELEX/jellyfish/MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6_realz_hamming1_signal_mer_counts_dumps.txt", col_names=FALSE) #224
kmer_background=read_delim("/scratch/project_2013895/SELEX/jellyfish/MEIS1_HT-SELEX_TGACCT20NGA_O_NTGACAN_1_6_realz_hamming1_background_mer_counts_dumps.txt", col_names=FALSE) #224

colnames(signal)=colnames(background)=colnames(kmer_signal)=colnames(kmer_background)=c("kmer","count")

test=signal %>% filter(kmer %in% kmer_signal$kmer)

ind=match(kmer_signal$kmer, signal$kmer)

test=cbind(kmer_signal,signal[ind,])
