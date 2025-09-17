#!/bin/bash


array=$1


outfolder=/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX_scaled_by_maxscore/
mkdir $outfolder

files=(/scratch/project_2013895/SELEX/streamed_Hamming12_kmers_exclude_Methyl-HT-SELEX/*.tsv)

echo ${#files[@]} #967


nro_files=${#files[@]} 


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_files-1)) ]] 
then
     echo $end_ind is greater than $(($nro_files-1))
     length=$(($nro_files-$start_ind+1))
fi


for file in "${files[@]:$start_ind:$length}"; do
  
  echo "$file"
  
  # Check that all lines have exactly 3 columns (tab-delimited)
  if ! awk -F'\t' 'NF != 3 { exit 1 }' "$file"; then
    echo "Skipping $file: not all lines have exactly 3 columns"
    continue
  fi
  
  
  cp "$file" $LOCAL_SCRATCH
  # Get max of column 2 (score)
  
  file_basename=$(basename "$file")
  max=$(awk -F'\t' 'BEGIN {max = 0} $2 > max { max = $2 } END { print max }' "$LOCAL_SCRATCH"/"$file_basename")



  # Output path
  out=$LOCAL_SCRATCH"/tmp_"$file_basename

  # Step 2: Scale column 2 by max
  awk -F'\t' -v OFS='\t' -v max="$max" '{
  norm = (max == 0) ? 0 : $2 / max
  print $1, norm, $3
  }' "$LOCAL_SCRATCH"/"$file_basename" > "$out"

  cp $out $outfolder$file_basename
  rm $LOCAL_SCRATCH"/"$file_basename
  rm $out
  
done

