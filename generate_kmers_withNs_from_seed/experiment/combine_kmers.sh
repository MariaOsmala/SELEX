#!/bin/bash

array=$1

outfolder=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/
mkdir $outfolder

realz_files=(/scratch/project_2013895/SELEX/kmers_from_fixed_N_seeds/*_canonical_realizations.txt)
hamming1_files=(/scratch/project_2013895/SELEX/kmers_from_fixed_N_seeds/*_canonical_hamming1.txt)
hamming2_files=(/scratch/project_2013895/SELEX/kmers_from_fixed_N_seeds/*_canonical_hamming2.txt)

echo ${#realz_files[@]} #3635
echo ${#hamming1_files[@]}
echo ${#hamming2_files[@]}

nro_files=${#hamming1_files[@]} 


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_files-1)) ]] 
then
     echo $end_ind is greater than $(($nro_files-1))
     length=$(($nro_files-$start_ind+1))
     end_ind=$(($nro_files-1))
fi

#start_ind=0
#end_ind=$(($nro_files-1))

for ((i=start_ind; i<=end_ind; i++)); do
  realz_file=${realz_files[i]}
  hamming1_file=${hamming1_files[i]}
  hamming2_file=${hamming2_files[i]}
  # do something with "$i" (the index) and "$file"
  #echo "$i -> $file"
  echo $i
  base=${hamming1_files[i]##*/}                # remove directories
  ids=${base%_canonical_hamming1.txt}       # remove fixed suffix
  
  #echo $realz_file
  #echo $hamming1_file
  echo $ids
  
  cat "$realz_file" "$hamming1_file" "$hamming2_file" > "${outfolder}${ids}.txt"
  cat "$realz_file" "$hamming1_file" > "${outfolder}${ids}_hamming1.txt"
  
  
done

