#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/seqkit/


module load seqkit
seqkit version


file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv" #3573
signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"
motifs="ID"

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


nro_pwms=${#signals[@]} #3573
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi

out_path=/scratch/project_2013895/SELEX/seqkit_kmer_counts_combined_fixed_N_seeds/

for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"
  # i=3179
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
  
  motif=${motifs[$i]}
  
  kmers=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/$motif".txt"
  kmers_fasta=/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/$motif".fasta"

  perl -ne 'print ">\n$_"' $kmers > $kmers_fasta #this is done already
 
  echo $input$sample_base
  
  echo $input$background_base
  
  
  rm $out_path$motif".txt"  
  touch $out_path$motif".txt"  
  while read kmer; do
    #echo $kmer
    count=$(seqkit locate -p "$kmer" --degenerate "$input$sample_base" 2>/dev/null | tail -n +2 | wc -l)
    echo -e "$kmer\t$count" >> $out_path$motif".txt"  
  done < $kmers
  
  kmer=NNNACGANNNNNNTCGTNNN
  count=$(seqkit locate -p "$kmer" --degenerate "$input$sample_base" 2>/dev/null | tail -n +2 | wc -l)
  # When searching by sequences, it's partly matching, and both positive and negative strands are searched.
  #-s, --by-seq search subseq on seq, both positive and negative strand are searched, and mismatch allowed using flag -m/--max-mismatch
  # -r, --use-regexp             patterns are regular expression 
  # -f, --pattern-file string    pattern file (one record per line)
  
  
  
  seqkit grep -s -r -f $kmers $input$sample_base
  
  seqkit grep -s -d -f $kmers $input$sample_base | seqkit seq -n | sort | uniq -c
  
  seqkit locate -f $kmers_fasta -d -j 8 $input$sample_base | tail -n +2 | cut -f2 | sort | uniq -c | awk '{print $2"\t"$1}' > counts.tsv


echo -e "sample\tpattern\tcount" > per_sample_counts.tsv
for fq in *.fastq.gz; do
  seqkit locate -f kmers.fa -d -j 8 "$fq" \
    | tail -n +2 | cut -f2 | sort | uniq -c \
    | awk -v s="$(basename "$fq")" '{print s"\t"$2"\t"$1}'
done >> per_sample_counts.tsv

seqkit locate --pattern-file $kmers_fasta --use-fmi -d -j 8 $input$sample_base > hits.tsv

#-s search in sequence, -i ignore case, -r regex, -v invert (keep reads without N).

zcat $sample    | awk 'END{print NR/4}'   # total
zcat $output$sample_base | awk 'END{print NR/4}'   # after filtering
zcat $sample | awk 'NR%4==2 && /[Nn]/ {c++} END{print c}'


done
