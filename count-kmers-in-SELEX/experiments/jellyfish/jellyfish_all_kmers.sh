#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/jellyfish

PATH=/projappl/project_2013895/softwares/jellyfish-2.3.1/bin:$PATH


file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv" #3573
signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"
motifs="ID"
seeds="seed"

filtered_files=/scratch/project_2013895/SELEX/data_exclude_reads_with_N/


readarray -t signals < <(
  awk -F'\t' -v c="$signal_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)
readarray -t backgrounds < <(
  awk -F'\t' -v c="$background_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

readarray -t motifs < <(
  awk -F'\t' -v c="$motifs" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

readarray -t seeds < <(
  awk -F'\t' -v c="$seeds" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)



#array=325 #MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2

nro_pwms=${#motifs[@]} #3573
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


jf_folder=/scratch/project_2013895/SELEX/jellyfish/jf_all_kmers/
mkdir -p $jf_folder
out_folder=/scratch/project_2013895/SELEX/jellyfish/filtered_counts_all_kmers/
mkdir -p $out_folder

if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi

for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"


  motif=${motifs[$i]}
  seed=${seeds[$i]}
  echo $seed
  
  seed_length=${#seed}
  
  background=${backgrounds[$i]}
  signal=${signals[$i]}
  
  signal_base="${signal##*/}"    
  signal_no_ext="${signal_base%.*}"
  
  background_base="${background##*/}"    
  background_no_ext="${background_base%.*}" 
  
  
  #filtered fastq files (reads with N removed)

  zcat $filtered_files$background_base | jellyfish count /dev/fd/0 -m $seed_length -s 100M -t 10 -C -o $jf_folder$motif"_background.jf" 
  zcat $filtered_files$signal_base | jellyfish count /dev/fd/0 -m $seed_length -s 100M -t 10 -C -o $jf_folder$motif"_signal.jf" 

  jellyfish dump -c $jf_folder$motif"_background.jf"  > $out_folder$motif"_background_counts.txt"
  jellyfish dump -c $jf_folder$motif"_signal.jf" > $out_folder$motif"_signal_counts.txt"
  
  kmers=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/$motif".txt"
  ./count_kmers_with_N.sh $kmers $out_folder$motif"_background_counts.txt" > counts.txt

done

