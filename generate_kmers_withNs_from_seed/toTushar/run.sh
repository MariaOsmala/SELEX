#Generate k-mers from seed
#Generates realizations, Hamming 1 neighbourhood, and Hamming 2 neighbourhood in separate files

#Remove possible old files
rm kmers/PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4*
rm kmer_counts/PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4*
#Arguments: Seed file_name output_folder
Rscript --no-save generate_kmers_with_fixed_Ns_hamming12.R NAAGRCGTCTTN PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4 kmers

# Combine the files
cat kmers/PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4_canonical_* > kmers/PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4_canonical_all.txt


out_path=kmer_counts/
mkdir -p $out_path
  
kmers=kmers/PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4_canonical_all.txt

#Before running this, I have removed reads containing Ns in SELEX signal and background data  
./count_kmers_canonical.sh -k $kmers -f SELEX_data/PROX1_eDBD_4_KR_TGCTCT40NTTA.fastq.gz -o $out_path"PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4_signal.txt" 
./count_kmers_canonical.sh -k $kmers -f SELEX_data/PROX1_eDBD_3_KR_TGCTCT40NTTA.fastq.gz  -o $out_path"PROX1_HT-SELEX_TGCTCT40NTTA_KR_NAAGRCGTCTTN_1_4_background.txt" 
 
 

