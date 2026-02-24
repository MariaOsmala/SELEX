#!/bin/bash
array=$1 #0-357

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/
#setwd("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments")
PATH=/projappl/project_2013895/softwares/spacek:$PATH
module load biokit 

#cd /scratch/project_2013895/SELEX/spacek/lambda_2026
#find . -type f -empty | wc -l #103
#remove empty files
#find . -type f -empty -delete #3470 left

file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/metadata_motifs_with_SELEX_data.tsv" #3573
signal_filenames="CSC_SELEX_filename"             
background_filenames="CSC_SELEX_background_filename"
motifs="ID"
seeds="seed"
ligands="ligand"


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

readarray -t ligands < <(
  awk -F'\t' -v c="$ligands" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $(h[c])}' "$file"
)

#array=325 #MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2

nro_pwms=${#motifs[@]} #3573
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     end_ind=$(($nro_pwms-1))
fi

#target="TFAP2C_MAX_CAP-SELEX_TTAGTC40NTCC_AY_TNSCCNNNGGSNNNNNNNNNNNNNNCACGTGN_1_3"
#for i in "${!motifs[@]}"; do [[ ${motifs[$i]} == "$target" ]] && echo $i; done
#for ((i=start_ind; i<=end_ind; i++)); do
for ((i=0; i<=3573; i++)); do
  #echo "i=$i"


  motif=${motifs[$i]}
  echo $motif
  seed=${seeds[$i]}
  ligand=${ligands[$i]}
  
  background=${backgrounds[$i]}
  sample=${signals[$i]}
  
  sample_base="${sample##*/}"    #MEIS1_TGACCT20NGA_O_6.fastq.gz
  sample_no_ext="${sample_base%.*}" #MEIS1_TGACCT20NGA_O_6.fastq
  
  background_base="${background##*/}"    
  background_no_ext="${background_base%.*}" 
  
  preprocessed_path=/scratch/project_2013895/SELEX/preprocessed_data/
  
  backgroud_test=$preprocessed_path"${background_no_ext%.fastq}".seq
  sample_test=$preprocessed_path"${sample_no_ext%.fastq}".seq
  
  out="${preprocessed_path}${background_no_ext}"
  [[ -f "$background_test" ]] || zcat -- "$background" > "$out"
  
  out="${preprocessed_path}${sample_no_ext}"
  [[ -f "$sample_test" ]] || zcat -- "$sample" > "$out"
  
  #wc -l ${preprocessed_path}${sample_no_ext} #320512
  #wc -l ${preprocessed_path}${background_no_ext} #252172
  
  
  out=$preprocessed_path"${sample_no_ext%.fastq}".fasta
  [[ -f "$sample_test" ]] || seqtk seq -A -- "${preprocessed_path}${sample_no_ext}" > "$out"
  out=$preprocessed_path"${background_no_ext%.fastq}".fasta
  [[ -f "$background_test" ]] || seqtk seq -A -- "${preprocessed_path}${background_no_ext}" > "$out"
  
  
  #wc -l $preprocessed_path"${sample_no_ext%.fastq}".fasta # 160256, numbers above divided by 2
  #wc -l $preprocessed_path"${background_no_ext%.fastq}".fasta # 126086
  
  [[ -f "$sample_test" ]] || seqkit stats -- $preprocessed_path"${sample_no_ext%.fastq}".fasta > /projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/seqkit_stats/${sample_no_ext%.fastq}.txt #80,128
  [[ -f "$background_test" ]] || seqkit stats -- $preprocessed_path"${background_no_ext%.fastq}".fasta > /projappl/project_2013895/SELEX/count-kmers-in-SELEX/Data/seqkit_stats/${background_no_ext%.fastq}.txt #63,043
  
  out=$preprocessed_path"${sample_no_ext%.fastq}".seq
  [[ -f "$sample_test" ]] || sed '/^>/d' -- $preprocessed_path"${sample_no_ext%.fastq}".fasta > "$out"
  
  out=$preprocessed_path"${background_no_ext%.fastq}".seq
  [[ -f "$backgroun_test" ]] || sed '/^>/d' -- $preprocessed_path"${background_no_ext%.fastq}".fasta > "$out"
  
  rm $preprocessed_path"${sample_no_ext%.fastq}".fasta
  rm $preprocessed_path"${background_no_ext%.fastq}".fasta
  rm ${preprocessed_path}${background_no_ext}
  rm ${preprocessed_path}${sample_no_ext}
  
  #Extract the sequence length from the ligand
  
  if [[ $ligand =~ ([0-9]+N) ]]; then
    ligand_length="${BASH_REMATCH[1]}"   
   
  fi
  
  
  out_folder=/scratch/project_2013895/SELEX/spacek/
  
  spacek40 "-$ligand_length" --f $preprocessed_path"${background_no_ext%.fastq}".seq $preprocessed_path"${sample_no_ext%.fastq}".seq $seed > $out_folder"lambda_2026/"$motif".txt"  #Produces also sample.seq_logo.svg (what info does this contain)
  
  svg=$preprocessed_path"${sample_no_ext%.fastq}".seq_logo.svg
  svg_path=/scratch/project_2013895/SELEX/spacek/svg/
  out=$svg_path"${motif}".svg
  #echo $svg
  [[ -f "$svg" ]] && mv -- $svg $out
  
  
done

#spacek40 -nogaps -c -20N background.seq sample.seq 10 10 1 > counts.txt
#spacek40 -c -20N background.seq sample.seq 10 10 1 > counts.txt


#Seed output
#            Total sequences	Hits	Fraction	Lower half 8mer total	Max 8mer count
#Background	63037	          1921	3	        742948	              424
#Signal    	80122	          11039	14	      838913	              3292
#Lambda 0.888

#So the lambda would be 
#(838913*63037)/(742948*80122)=0.8883871



# Generate SVG logo file from input PWM, also generates .png if convert is installed in path
# module load imagemagick
#./spacek40 --logo [PWM file name] [output file name (opt)]

# Barcode logo generation usage:
#./spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname <pfm_file_name>
#./spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm

# Refer TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm as a pfm_file format example.
#	Row 1: A
#	Row 2: C
#	Row 3: G
#	Row 4: T

#Example output file: TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm.svg

#Generate difference logo (PWM1-PWM2)
#./spacek40 --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] 

#Calculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. 
#Default kmers are ungapped, default lengths are PWM widths (8 for gapped)
#./spacek40 --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] 

# Options:
# 
# -paths	Generate PWM logos using paths and not fonts
# -neg	Allow negative values in PWM (using paths)
# -nocall	 do not include kmer-based prediction calls into logos
# -barcodelogo	Generate PWM barcode logos using colored rectangles
# -heightscaledbars	Scale height of barcode logos based on nucleotide frequency
# -maxheightscaledbars	Scale height of barcode logos based on maximum nucleotide frequency at position
# -colorscaledbars	Scale color of barcode logos based on nucleotide frequency
# -label	Include white IUPAC label to barcodelogo (cutoff = more than 50% of max)
# -neg	Allow negative values in PWM (using paths)
# -noname	Do not include filename in the logo
# -editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)	print table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance
# -longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]
# -iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]
# -match=[kmer1],[kmer2]	Include spacing and orientation heatmap for these kmers in kmer summary svg
# -hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]
# -q	print frequencies instead of counts
# -n	print nucleotide counts for full sequences
# -s	Print values for different gap lengths on same line
# -c	Print raw counts from both files
# -u	Use only unique input sequences
# -p	Print p-values (Winflat program needs to be installed in path)
# -i	Print incidence of all input sequences
# -e	Use even background instead of background from file
# -lim=[position followed by strand (one or both must be given)]	 Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)
# -exc=[position][strand],[position][strand]‚..,[position][strand]	Exclude indicated positions and strands (e.g. -exc=4,3F,6R)
# -m=[number]	Generate PWM using multinomial [number] distribution
# -dimer	print counts for dimeric sequences only
# -kmer	Print kmer count table
# -dinuc	print dinucleotide data to output and svg
# -ic	Output information content for all spacings
# -nogaps	Count only kmers without gaps
# -allgaps	Count kmers with gaps in any position (default is only middle 1 or 2 positions)
# -mono	Use mononucleotide background instead of multinomial background
# -both	Count both instances of palindromic hits
# -bothifnotequal	Count both instances of palindromic hits if hit scores are not equal
# -forwardonly	Count only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)
# -reverseonly	Count only reverse instance of palindromic hits
# -forwardifequal	Count only forward instance of palindromic hits if scores equal (default alternates between strands)
# -reverseifequal	Count only reverse instance of palindromic hits if scores equal
# -14N	Set sequence length to 14N (default is 20N)
# -30N	Set sequence length to 30N
# -40N	Set sequence length to 40N
# -fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)
# -mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits
# -kl=[kmer length] kmer length used to calculate lambda (default 8)
# -x	Extended output for debugging
# -o=[output file name] name of svg output file



# Autoseed usage:
# $ ./totalautoseed -40N <Background sequence> <Signal sequence> <huddinge_distance> 
#  <min_seed_length> <max_seed_length> <length_cutoff> <iupac_cutoff> <local_max_cutoff_count> <logo_number_cutoff>
 
# Ex) In case 20N SELEX oligos.
# $ ./totalautoseed -20N TTCTTC20NTTCKH3_sig.seq TTCTTC20NTTCKH4_sig.seq 1 8 10 0.35 - 50 100
# TTCTTC20NTTCKH3_sig.seq and TTCTTC20NTTCKH4_sig.seq are 20 base SELEX sequence files (without headers).
 
# Ex) In case 40N SELEX oligos.
# $ ./totalautoseed -40N aaa.txt bbb.txt 1 8 10 0.35 - 50 100
# aaa.txt and bbb.txt are 40 base SELEX sequence files (without headers).
 
# Example output file:
# TTCTTC20NTTCKH4_sig.seq_1_8_10_0.35_-_50_logos.svg


#---
#Dinucleotide heatmap generation usage:
#$ ./spacek40 --f -dinuc -nocall -m=<multinomial> <Nlength_option> <Background sequence> <Signal sequence> <Seed> <cutoff>

#Ex) In case 20N SELEX oligos.
#$ ./spacek40 --f -dinuc -nocall -m=1 -20N TTCTTC20NTTCKH3_sig.seq TTCTTC20NTTCKH4_sig.seq NATGACGTSATN 200
#TTCTTC20NTTCKH3_sig.seq and TTCTTC20NTTCKH4_sig.seq are 20 base SELEX sequence files (without headers).

#Ex) In case 40N SELEX oligos.
#$ ./spacek40 --f -dinuc -nocall -m=1 -40N aaa.txt bbb.txt NCCGGAAN 200
#aaa.txt and bbb.txt are 40 base SELEX sequence files (without headers).

#Example output file:
#TTCTTC20NTTCKH4_sig.seq_logo.svg


