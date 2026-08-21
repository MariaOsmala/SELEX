#!/bin/bash
array=$1 #0-455

#218-> requires more memory than 50G

cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments/
PATH=/projappl/project_2013895/softwares/jellyfish-2.3.1/bin:$PATH

module load biokit
file="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/unique_SELEX.txt" #4552


mapfile -t files < $file


nro_files=${#files[@]} #4552
start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_files-1)) ]] 
then
     echo $end_ind is greater than $(($nro_files-1))
     end_ind=$(($nro_files-1))
fi




filtered_files=/scratch/project_2013895/SELEX/data_exclude_reads_with_N/

#jf_folder=/scratch/project_2013895/SELEX/lambda_computation/jf/
jf_folder=/scratch/project_2013895/SELEX/lambda_computation/jf_without_C/

#out_folder=/scratch/project_2013895/SELEX/lambda_computation/filtered_counts/
out_folder=/scratch/project_2013895/SELEX/lambda_computation/filtered_counts_without_C/

mkdir -p $jf_folder
mkdir -p $out_folder

i=$(printf '%s\n' "${files[@]}" | grep -n 'FOXO1_TTCAGC20NTA_AF_2' | cut -d: -f1)
# grep -n gives 1-based line numbers, so subtract 1 for bash array index:
i=$((i - 1))
echo "$i"



for ((i=start_ind; i<=end_ind; i++)); do
  echo "i=$i"
  #i=2280

  file=${files[$i]}
  base="${file##*/}"    
  base_no_ext="${base%.fastq.gz*}" 

  #filtered fastq files (reads with N removed)
  # Step 1: Count 8-mers with jellyfish
  
  #The -C flag (canonical/both strands) — do you want to count canonical k-mers (collapsing reverse complements), 
  #or count each strand separately? 
  
  #For SELEX data you may want both strands counted independently (omit -C).
  
  rm $jf_folder$base_no_ext".jf"
  rm $out_folder$base_no_ext".txt"
  
  n_reads_orig=$(( $(zcat $file | wc -l) / 4 ))
  echo "orig reads: $n_reads_orig"
  
  n_reads_after_N_filtering=$(( $(zcat $filtered_files$base | wc -l) / 4 ))
  echo "reads after N filtering: $n_reads_after_N_filtering"
  
  # --- 2. Produce reverse-complemented copy ---
  seqkit seq -r -p -t DNA $filtered_files$base | gzip > $filtered_files"rev_comp_"$base
  

  jellyfish count -m 8 -s 100M -t 8 -o $jf_folder$base_no_ext".jf"  \
  <(pigz -dc $filtered_files$base) <(pigz -dc $filtered_files"rev_comp_"$base)


  jellyfish dump -c $jf_folder$base_no_ext".jf" > $out_folder$base_no_ext".txt"
  


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


