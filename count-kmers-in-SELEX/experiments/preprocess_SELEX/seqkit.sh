#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/


module load biokit


file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv" #3573
signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"

output=/scratch/project_2013895/SELEX/data_exclude_reads_with_N/

readarray -t signals < <(
  awk -F'\t' -v c="$signal_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)
readarray -t backgrounds < <(
  awk -F'\t' -v c="$background_filenames" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)



nro_pwms=${#signals[@]} #3573
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi

for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"
  
  background=${backgrounds[$i]}
  sample=${signals[$i]}
  
  #sample=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/BARHL2_TCCAGT40NGAC_AI_3.fastq.gz
  #background=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/BARHL2_TCCAGT40NGAC_AI_2.fastq.gz
  
  sample_base="${sample##*/}"    
  sample_no_ext="${sample_base%.*}"
  
  background_base="${background##*/}"    
  background_no_ext="${background_base%.*}" 
  backgroud_test=$output$background_base
  sample_test=$output$sample_base
  
 
  [[ -f "$sample_test" ]] || seqkit grep -s -i -r -v -p N -- $sample | gzip -c > $output$sample_base
  
  [[ -f "$background_test" ]] || seqkit grep -s -i -r -v -p N -- $background | gzip -c > $output$background_base



#-s search in sequence, -i ignore case, -r regex, -v invert (keep reads without N).

zcat $sample    | awk 'END{print NR/4}'   # total
zcat $output$sample_base | awk 'END{print NR/4}'   # after filtering
zcat $sample | awk 'NR%4==2 && /[Nn]/ {c++} END{print c}'


done
