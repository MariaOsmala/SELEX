#!/bin/bash


array=$1

outfolder=/scratch/project_2013895/SELEX/streamed_kmers_scaled/
mkdir $outfolder

files=(/scratch/project_2013895/SELEX/streamed_kmers/*.tsv)

echo ${#files[@]} #3774 -18(wrong format), there was 202 failed experiments


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
  # Get min and max of column 2 (score)
  
  file_basename=$(basename "$file")
  read min max <<< $(awk 'NR==1 { min = max = $2 } $2 < min { min = $2 } $2 > max { max = $2 } END { print min, max }' "$LOCAL_SCRATCH"/"$file_basename")

  # Output path
  out=$LOCAL_SCRATCH"/tmp_"$file_basename


  # Scale and write to output
  awk -v min="$min" -v max="$max" -F'\t' -v OFS='\t' '{
    scaled = (max == min) ? 1 : ($2 - min) / (max - min)
    print $1, scaled, $3
  }' "$LOCAL_SCRATCH"/"$file_basename" > "$out"
  
  cp $out $outfolder$file_basename
  rm $LOCAL_SCRATCH"/"$file_basename
  rm $out
  
done

