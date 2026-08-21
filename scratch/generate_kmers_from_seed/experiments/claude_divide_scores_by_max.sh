#!/bin/bash


array=$1


outfolder=/scratch/project_2013895/SELEX/scored_scaled_kmers_from_long_degenerate_seeds/
mkdir $outfolder


files=(/scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/*.tsv)

echo ${#files[@]} #9782 (0-97)


nro_files=${#files[@]} 


start_ind=$(($array*100)) #100
end_ind=$((($array+1)*100 -1)) #100
length=100 #100


if [[ $end_ind -gt $(($nro_files-1)) ]] 
then
     echo $end_ind is greater than $(($nro_files-1))
     end_ind=$(($nro_files-1))
fi


for file in "${files[@]:$start_ind:$length}"; do
  
  echo "$file"
  
  #awk -F'\t' 'NR==1{print NF; exit}' /scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3.tsv #3
  #awk -F'\t' 'NR==1{print NF; exit}' /scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3_SELEX_counts.tsv #11
  #awk -F'\t' 'NR==1{print NF; exit}' /scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3_jellyfish_counts.tsv #7
  #awk -F'\t' 'NR==1{print NF; exit}' /scratch/project_2013895/SELEX/scored_kmers_from_long_degenerate_seeds/Alx1_HT-SELEX_TAAAGC20NCG_Z_NNYAATTANN_1_3_combined_counts.tsv #15
  
  ncols=$(awk -F'\t' '$0 !~ /^#/ && NF { print NF; exit }' "$file")
  
  cp "$file" $LOCAL_SCRATCH
  # Get max of column 2 (score)
  
  file_basename=$(basename "$file")
  
  case "$ncols" in
  3)  echo "Do thing for 3 columns"
  
  in="$LOCAL_SCRATCH/$file_basename"
  max=$(awk -F'\t' 'BEGIN {max = 0} $2 > max { max = $2 } END { print max }' "$in") 
  # Output path
  out=$LOCAL_SCRATCH"/tmp_"$file_basename
  
  awk -F'\t' -v OFS='\t' -v max="$max" '
  {
    norm = (max == 0 || $2 == "") ? "" : ($2 / max)
    # keep col2, then add scaled value, then the rest
    print $1, $2, norm, $3
  }' "$in" > "$out"

  cp $out $outfolder$file_basename
  rm $LOCAL_SCRATCH"/"$file_basename
  rm $out ;; #.tsv
  7)  echo "Do thing for 7 columns (jellyfish)"  
      IFS=$'\t' read -r -a COLS < <(head -n1 $file)
      printf '%s\n' "${COLS[@]}"
      
      #need to scale motif_match_score and jellyfish Corrected count, add the normalized scores
      
      in=$LOCAL_SCRATCH"/"$file_basename
      out=$LOCAL_SCRATCH"/tmp_"$file_basename
      
      awk -F'\t' -v OFS='\t' \
        -v COLA='motif_match_score' \
        -v COLB='jellyfish Corrected count' \
        -v LABA='motif_match_score_scaled' \
        -v LABB='jellyfish Corrected count scaled' '
      # Pass 1: locate columns, gather maxima
      FNR==1 {
        m=j=0
        for (i=1;i<=NF;i++) {
          h=$i; gsub(/^"+|"+$/,"",h)
          if (h==COLA) m=i
          if (h==COLB) j=i
        }
        if (NR==1) next            # skip header only on first file (pass 1)
      }
      NR==FNR {
        if (m && $m!="") {v=$m+0; if (v>maxm) maxm=v}
        if (j && $j!="") {v=$j+0; if (v>maxj) maxj=v}
        next
      }
      
      # Pass 2: print header with inserted labels, then rows with inserted scaled values
      FNR==1 {
        if (!m || !j) { print "ERROR: required columns not found" > "/dev/stderr"; exit 1 }
        for (i=1;i<=NF;i++) {
          printf "%s", $i
          if (i==m) printf "%s%s", OFS, LABA
          if (i==j) printf "%s%s", OFS, LABB
          printf (i<NF?OFS:ORS)
        }
        next
      }
      {
        a = (m && maxm>0 && $m!="") ? ($m+0)/maxm : ""
        b = (j && maxj>0 && $j!="") ? ($j+0)/maxj : ""
        # Uncomment for fixed decimals:
        # if (a!="") a=sprintf("%.6f", a); if (b!="") b=sprintf("%.6f", b)
      
        for (i=1;i<=NF;i++) {
          printf "%s", $i
          if (i==m) printf "%s%s", OFS, a
          if (i==j) printf "%s%s", OFS, b
          printf (i<NF?OFS:ORS)
        }
      }
      ' "$file" "$file" > "$out"
     cp $out $outfolder$file_basename
     rm $LOCAL_SCRATCH"/"$file_basename
     rm $out ;; # jellyfish
  11) echo "Do thing for 11 columns (SELEX)" 
      IFS=$'\t' read -r -a COLS < <(head -n1 $file)
      printf '%s\n' "${COLS[@]}"
      
      #need to scale motif_match_score and jellyfish Corrected count, add the normalized scores
      
      in=$LOCAL_SCRATCH"/"$file_basename
      out=$LOCAL_SCRATCH"/tmp_"$file_basename
      
      awk -F'\t' -v OFS='\t' \
        -v COLA='motif_match_score' \
        -v COLB='Corrected count' \
        -v LABA='motif_match_score_scaled' \
        -v LABB='Corrected count scaled' '
      # Pass 1: locate columns, gather maxima
      FNR==1 {
        m=j=0
        for (i=1;i<=NF;i++) {
          h=$i; gsub(/^"+|"+$/,"",h)
          if (h==COLA) m=i
          if (h==COLB) j=i
        }
        if (NR==1) next            # skip header only on first file (pass 1)
      }
      NR==FNR {
        if (m && $m!="") {v=$m+0; if (v>maxm) maxm=v}
        if (j && $j!="") {v=$j+0; if (v>maxj) maxj=v}
        next
      }
      
      # Pass 2: print header with inserted labels, then rows with inserted scaled values
      FNR==1 {
        if (!m || !j) { print "ERROR: required columns not found" > "/dev/stderr"; exit 1 }
        for (i=1;i<=NF;i++) {
          printf "%s", $i
          if (i==m) printf "%s%s", OFS, LABA
          if (i==j) printf "%s%s", OFS, LABB
          printf (i<NF?OFS:ORS)
        }
        next
      }
      {
        a = (m && maxm>0 && $m!="") ? ($m+0)/maxm : ""
        b = (j && maxj>0 && $j!="") ? ($j+0)/maxj : ""
        # Uncomment for fixed decimals:
        # if (a!="") a=sprintf("%.6f", a); if (b!="") b=sprintf("%.6f", b)
      
        for (i=1;i<=NF;i++) {
          printf "%s", $i
          if (i==m) printf "%s%s", OFS, a
          if (i==j) printf "%s%s", OFS, b
          printf (i<NF?OFS:ORS)
        }
      }
      ' "$file" "$file" > "$out"
     cp $out $outfolder$file_basename
     rm $LOCAL_SCRATCH"/"$file_basename
     rm $out ;; #SELEX_counts
  15) echo "Do thing for 15 columns" 
  IFS=$'\t' read -r -a COLS < <(head -n1 $file)
      printf '%s\n' "${COLS[@]}"
      
      #need to scale motif_match_score and Corrected count, add the normalized scores
      
      in=$LOCAL_SCRATCH"/"$file_basename
      out=$LOCAL_SCRATCH"/tmp_"$file_basename
      
      awk -F'\t' -v OFS='\t' \
      -v COLA='motif_match_score' \
      -v COLB='Corrected count' \
      -v COLC='jellyfish Corrected count' \
      -v LABA='motif_match_score_scaled' \
      -v LABB='Corrected count scaled' \
      -v LABC='jellyfish Corrected count scaled' '
    # ---------- PASS 1: locate columns & collect maxima ----------
    FNR==1 {
      if (NR==1) {
        m=j=k=0
        for (i=1;i<=NF;i++) {
          h=$i; gsub(/^"+|"+$/,"",h)       # strip optional quotes
          if (h==COLA) m=i
          if (h==COLB) j=i
          if (h==COLC) k=i
        }
        next
      }
    }
    NR==FNR {
      if (m && $m!="") {v=$m+0; if (v>maxm) maxm=v}
      if (j && $j!="") {v=$j+0; if (v>maxj) maxj=v}
      if (k && $k!="") {v=$k+0; if (v>maxk) maxk=v}
      next
    }
    
    # ---------- PASS 2: print header + data with inserted scaled ----------
    FNR==1 {
      if (!m || !j || !k) {
        print "ERROR: one of the required columns not found" > "/dev/stderr"; exit 1
      }
      for (i=1;i<=NF;i++) {
        printf "%s", $i
        if (i==m) printf "%s%s", OFS, (LABA!=""?LABA:COLA" scaled")
        if (i==j) printf "%s%s", OFS, (LABB!=""?LABB:COLB" scaled")
        if (i==k) printf "%s%s", OFS, (LABC!=""?LABC:COLC" scaled")
        printf (i<NF?OFS:ORS)
      }
      next
    }
    {
      a = (maxm>0 && $m!="") ? (($m+0)/maxm) : ""
      b = (maxj>0 && $j!="") ? (($j+0)/maxj) : ""
      c = (maxk>0 && $k!="") ? (($k+0)/maxk) : ""
      # For fixed decimals: if(a!="")a=sprintf("%.6f",a); if(b!="")b=sprintf("%.6f",b); if(c!="")c=sprintf("%.6f",c)
    
      for (i=1;i<=NF;i++) {
        printf "%s", $i
        if (i==m) printf "%s%s", OFS, a
        if (i==j) printf "%s%s", OFS, b
        if (i==k) printf "%s%s", OFS, c
        printf (i<NF?OFS:ORS)
      }
    }
    ' "$in" "$in" > "$out"
    cp $out $outfolder$file_basename
    rm $LOCAL_SCRATCH"/"$file_basename
    rm $out ;; #combined_counts
  *)  echo "Unexpected column count: $ncols" ;;
esac




  
done

