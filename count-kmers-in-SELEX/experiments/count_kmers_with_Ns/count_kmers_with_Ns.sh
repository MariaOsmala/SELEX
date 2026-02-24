#!/bin/bash


array=$1
addition=$2

file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv" #3573
signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"
motifs="ID"

array=$(($array+$addition)) #100

input=/scratch/project_2013895/SELEX/data_exclude_reads_with_N/


readarray -t signals < <(
  awk -F'\t' -v c="$signal_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)
readarray -t backgrounds < <(
  awk -F'\t' -v c="$background_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

readarray -t motifs < <(
  awk -F'\t' -v c="$motifs" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

#1500-1509,2270-2279

nro_pwms=${#signals[@]} #3573
#start_ind=$(($array*10)) #100
#end_ind=$((($array+1)*10 -1)) #100

start_ind=$(($array*1)) #100
end_ind=$((($array+1)*1 -1)) #100



if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi

out_path=/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/

for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"
  #i=3179
  #i=2280
  background=${backgrounds[$i]}
  sample=${signals[$i]}
  
  sample_base="${sample##*/}"    
  sample_no_ext="${sample_base%.*}"
  sample_no_ext="${sample_no_ext%.*}"
  
  background_base="${background##*/}"    
  background_no_ext="${background_base%.*}" 
  
  motif=${motifs[$i]}
  
  kmers=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/$motif".txt"
  #kmers=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/$motif"_hamming1.txt"
  
  #chmod +x count_kmers.sh
  # tested for HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_NNNACGANNNNNNTCGTNNN_1_3u
  #time ./count_kmers.sh -k $kmers -f $input$sample_base -o $out_path$motif".txt" #  #this does not count also reverse complements 
  #real    0m59.752s 
  #user    0m58.006s
  #sys     0m1.412s
  #time ./count_kmers_faster.sh $kmers $input$sample_base > $out_path$motif"_singlepass.txt"
  #real    13m0.524s
  #user    12m58.516s
  #sys     0m0.182s
  #time ./count_kmers_canonical.sh -k $kmers -f $input$sample_base -o $out_path$motif"_canonical_hamming1.txt"
  #./count_kmers_canonical.sh -k $kmers -f $input$sample_base -o $out_path$motif"_hamming1_signal.txt" 
  ./count_kmers_canonical.sh -k $kmers -f $input$sample_base -o $out_path$motif"_signal.txt" 
  ./count_kmers_canonical.sh -k $kmers -f $input$background_base -o $out_path$motif"_background.txt" 
  #real    1m5.318s
  #user    1m3.458s
  #sys     0m1.625s
 done 
 
 