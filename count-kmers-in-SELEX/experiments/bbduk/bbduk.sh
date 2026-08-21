#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/bbduk

module load biojava
/projappl/project_2013895/softwares/bbmap


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

out_path=/scratch/project_2013895/SELEX/bbduk_kmer_counts_combined_fixed_N_seeds/

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
 
 #Compares reads to the kmers in a reference dataset, optionally 
 #allowing an edit distance. Splits the reads into two outputs - those that 
 #match the reference, and those that don't. Can also trim (remove) the matching 
 #parts of the reads rather than binning the reads.
 #Please read bbmap/docs/guides/BBDukGuide.txt for more information.

 #Usage:  bbduk.sh in=<input file> out=<output file> ref=<contaminant files>

 #Input may be stdin or a fasta or fastq file, compressed or uncompressed.
 #If you pipe via stdin/stdout, please include the file type; e.g. for gzipped 
 #fasta input, set in=stdin.fa.gz

  ./bbduk.sh in=$input$sample_base ref=$kmers_fasta k=20 mm=f rcomp=t out=$output$motif".fastq"


done
