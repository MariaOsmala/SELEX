#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/


module load biokit

file="unique_SELEX.txt" #4552
output=/scratch/project_2013895/SELEX/data_exclude_reads_with_N/

mapfile -t files < $file


nro_files=${#files[@]} #3573
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_files-1)) ]] 
then
     echo $end_ind is greater than $(($nro_files-1))
     end_ind=$(($nro_files-1))
fi


for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"
  echo ${files[$i]}
  
  file=${files[$i]}
  base="${file##*/}"    
 #Take $file, remove all sequence records whose sequence contains N (or n)
 # and write the remaining records as a gzipped file to $output$base.
 #One subtle point: with -r, pattern N is just the regex N, 
 #so here it behaves the same as a plain literal N. If the input is FASTQ, entire reads are removed
  seqkit grep -s -i -r -v -p N -- $file | gzip -c > $output$base
  
  #-s search in sequence, -i ignore case, -r regex, -v invert (keep reads without N).

  zcat $file  | awk 'END{print NR/4}'   # total
  zcat $output$base | awk 'END{print NR/4}'   # after filtering
  zcat $file | awk 'NR%4==2 && /[Nn]/ {c++} END{print c}' #counts how many FASTQ sequence lines contain N or n.


done
