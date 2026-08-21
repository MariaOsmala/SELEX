#!/bin/bash
array=$1 #0-267

module load r-env/452

output_tables="/scratch/project_2013895/SELEX/lambda_computation/output_tables/"
output="/scratch/project_2013895/SELEX/lambda_computation/output/"

file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments/unique_signal_and_background.tsv"

signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"


readarray -t signals < <(
  awk -F'\t' -v c="$signal_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)
readarray -t backgrounds < <(
  awk -F'\t' -v c="$background_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

readarray -t ligands < <(
  awk -F'\t' -v c="ligand" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

nro_pwms=${#signals[@]} #2681
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


 # i=$(printf '%s\n' "${signals[@]}" | grep -n 'IRX5_TBX20_3_YWIIII_TTGGCG40NGGGA' | cut -d: -f1)
 # # grep -n gives 1-based line numbers, so subtract 1 for bash array index:
 # i=$((i - 1))
 # echo "$i"


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi


#for ((i=0; i<=2680; i++)); do
for ((i=start_ind; i<=end_ind; i++)); do

  echo "i=$i"

  ligand=${ligands[$i]}
  [[ $ligand =~ ([0-9]+) ]] && ligand_length="${BASH_REMATCH[1]}"
  
  background=${backgrounds[$i]}
  signal=${signals[$i]}
  
  signal_base="${signal##*/}"    
  signal_no_ext="${signal_base%".fastq.gz"*}" 
  
  background_base="${background##*/}"    #fastq.gz
  background_no_ext="${background_base%".fastq.gz"*}" 
  
  Rscript ../code/compute_lambda_R.R --sig $signal --bg  $background --length $ligand_length --kmer 8 --out $output_tables$signal_no_ext"_"$background_no_ext".tsv" > $output$signal_no_ext"_"$background_no_ext".tsv"
  
  
done