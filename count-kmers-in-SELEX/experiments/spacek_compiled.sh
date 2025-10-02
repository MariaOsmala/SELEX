
cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments
#setwd("/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments")
spacek_path=/projappl/project_2013895/softwares/spacek

PATH=/projappl/project_2013895/softwares/spacek:$PATH

#Usage: ./spacek --tool -options [background file name] [sample file name] 
#        [shortest kmer] [longest kmer] [minimum count or fold (with option -s) for printing] [PWM file (only with option -kmer)]

#Count number of occurrences of subsequences (kmers) of the indicated size. 
#Indicates if count is higher than sequences within 1 Huddinge distance.
#Writes also a graphical summary to an svg file.         
        

Other tools:
--f [background file name or - for none] [sample file name] [seed kmer sequence (consensus or IUPAC) or 
     number for automatic (number specifies rank of local max start seed)] [minimum incidence] or
--pwm [background file name or - for none] [sample file name] [pwm file name] [score cut-off e.g. 0.1 or 10%] [minimum incidence] 
        Find occurrences of input sequence / pwm and print flanking preference, spacing and orientation information. 
        Graphical summary of data is also written to an svg file
 --logo [PWM file name] [output file name (opt)]
        Generate SVG logo file from input PWM, also generates .png if convert is installed in path
  --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] 
        Generate difference logo (PWM1-PWM2)
  --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] 
        Calculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. 
        Default kmers are ungapped, default lengths are PWM widths (8 for gapped)
  --pwmalign [PWM1 file name] [PWM2 file name] [PWM3 file name (optional)]
        Align PWMs and generate a svg with logos aligned to PWM2, with boxes indicating similar base positions.
        Option -contacts shows also base contacts of PWM1 and PWM3 on PWM2.


Options:

Logo and sequence formatting:
  -contacts     Show hydrogen bond contacts (PWMs where contact information is included must be used)
  -paths        Generate PWM logos using paths and not fonts
  -neg  Allow negative values in PWM (using paths)
  -nocall        do not include kmer-based prediction calls into logos
  -noname       Do not include filename in the logo
  -barcodelogo  Generate PWM barcode logos using colored rectangles
  -heightscaledbars     Scale height of barcode logos based on nucleotide frequency
  -maxheightscaledbars  Scale height of barcode logos based on maximum nucleotide frequency at position
  -colorscaledbars      Scale color of barcode logos based on nucleotide frequency
  -label        Include white IUPAC label to barcodelogo (cutoff = more than 50% of max)
  -rna  Input sequence is from RNA, use uracil (U) instead of thymidine (T) in sequences and logos
  -CpG  Highlight CpGs in logos if they occur with lower (opaque) or higher (black outline) frequency in sample file than in control file 

Input sequences:
  -14N  Set sequence length to 14N (default is 20N)
  -30N  Set sequence length to 30N
  -40N  Set sequence length to 40N
  -u    Use only unique input sequences

Kmer counting:
  -nogaps       Count only kmers without gaps
  -allgaps      Count kmers with gaps in any position (default is only middle 1 or 2 positions)
  -longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]
  -iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]

Matching and PWM generation:
  -m=[number]   Generate PWM using multinomial [number] distribution
  -iterate      iterate input seed based on generated one hit pwm allowing up to 2 bp longer seed per round
  -iterate-samesize     iterate seed, limit seed size to specified length
  -iterate-fast  iterate seed using free seed length
  -lim=[position followed by strand (one or both must be given)]         Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)
  -exc=[position][strand],[position][strand]‚..,[position][strand]      Exclude indicated positions and strands (e.g. -exc=4,3F,6R)
  -both Count both instances of palindromic hits
  -bothifnotequal       Count both instances of palindromic hits if hit scores are not equal
  -forwardonly  Count only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)
  -reverseonly  Count only reverse instance of palindromic hits
  -forwardifequal       Count only forward instance of palindromic hits if scores equal (default alternates between strands)
  -reverseifequal       Count only reverse instance of palindromic hits if scores equal
  -e    Use even background instead of background from file
  -mono Use mononucleotide background instead of multinomial background
  -fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)
  -mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits
  -kl=[kmer length] kmer length used to calculate lambda and automatic start seed (default 8)
  
Text output:
  -q    Print frequencies instead of counts
  -n    Print nucleotide counts for full sequences
  -s    Print values for different gap lengths on same line
  -c    Print raw counts from both input files
  -p    Print p-values (Winflat program needs to be installed in path)
  -i    Print incidence of all input sequences
  -kmer Print kmer count table with scores for each kmer against the PWM indicated
  -editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)      Print table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance
  -dimer        Print counts for dimeric sequences only
  -ic   Output information content for all spacings
  -dinuc        Print dinucleotide data to output and svg (counting dinucleotides uses multinomial 2)
  -x    Extended output for debugging

Svg output:
  -o=[output file name] base name of svg output file(s)
  -match=[kmer1],[kmer2]        Include spacing and orientation heatmap for these kmers in kmer summary svg
  -hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]
  -xyplot       Generate scatterplot of kmer counts for the shortest kmer, with CpG-containing kmers indicated in red
  -eoplot       Generate scatterplot of kmer counts observed and expected from the one hit PWM
  
  



#k-mers from seed here scored
#/scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore

motif=HOXA10_HT-SELEX_TAGGAT30NAGT_AI_NGTCGTWAAANN_1_4 
seed=NGTCGTWAAANN
sample=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/HOXA10_TAGGAT30NAGT_AI_4.fastq.gz	
background=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/HOXA10_TAGGAT30NAGT_AI_3.fastq.gz


motif=MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2
seed=NNCACGTGNN
kmer_sequence=ACCACGTGCT
background=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/MAX_TGACCT20NGA_Y_1.fastq.gz 
sample=/scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp/MAX_TGACCT20NGA_Y_2.fastq.gz 


module load biokit 
zcat $background > background.fq
zcat $sample > sample.fq

wc -l sample.fq # 320512
wc -l background.fq # 252172



seqtk seq -A background.fq > background.fasta
seqtk seq -A sample.fq > sample.fasta

wc -l sample.fasta # 160256
wc -l background.fasta # 126086


seqkit stats sample.fasta #80,128
seqkit stats background.fasta #63,043


# file              format  type  num_seqs     sum_len  min_len  avg_len  max_len
# background.fasta  FASTA   DNA    596,116  17,883,480       30       30       30
# file          format  type  num_seqs     sum_len  min_len  avg_len  max_len
# sample.fasta  FASTA   DNA    459,522  13,785,660       30       30       30


sed '/^>/d' sample.fasta > sample.seq
sed '/^>/d' background.fasta > background.seq


#spacek40  -30N --f background.seq sample.seq $seed $mi > seed_output.txt
#spacek40  -30N --f background.seq sample.seq $kmer_sequence $mi > kmer_output.txt 

spacek40  -20N --f background.seq sample.seq $seed > seed_output.txt #Produces also sample.seq_logo.svg (what info does this contain)

spacek40 -nogaps background.seq sample.seq 10 10 1 > counts.txt #752465 


spacek40 -nogaps -c -20N background.seq sample.seq 10 10 1 > counts.txt
spacek40 -c -20N background.seq sample.seq 10 10 1 > counts.txt


#Seed output
            Total sequences	Hits	Fraction	Lower half 8mer total	Max 8mer count
Background	63037	          1921	3	        742948	              424
Signal    	80122	          11039	14	      838913	              3292
Lambda 0.888

#So the lambda would be 
(838913*63037)/(742948*80122)=0.8883871


EOF encountered in file 0 on line 63038
Number of sequences in file 0: 63037

EOF encountered in file 1 on line 80123
Number of sequences in file 1: 80122



EOF encountered in file 0 on line 510563
Number of sequences in file 0: 510562

EOF encountered in file 1 on line 365858
Number of sequences in file 1: 365857


sed -n '510562,510564p' background.seq

sed -n '365857,365859p' sample.seq

tr -d '\032' < background.seq > clean_background.seq
tr -d '\032' < sample.seq > clean_sample.seq

tr -cd '\11\12\15\40-\176' < clean_background.seq > background.seq
tr -cd '\11\12\15\40-\176' < clean_sample.seq > sample.seq

spacek40  -30N --f clean_background.seq clean_sample.seq $seed $mi > seed_output.txt

# Generate SVG logo file from input PWM, also generates .png if convert is installed in path
# module load imagemagick
./spacek40 --logo [PWM file name] [output file name (opt)]

# Barcode logo generation usage:
./spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname <pfm_file_name>
./spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm

# Refer TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm as a pfm_file format example.
#	Row 1: A
#	Row 2: C
#	Row 3: G
#	Row 4: T

#Example output file: TAAAAT20NCG_KD_ARAMCACGTGKTYT_m2_c4_short.pfm.svg

#Generate difference logo (PWM1-PWM2)
./spacek40 --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] 

#Calculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. 
#Default kmers are ungapped, default lengths are PWM widths (8 for gapped)
./spacek40 --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] 


Options:

-paths	Generate PWM logos using paths and not fonts
-neg	Allow negative values in PWM (using paths)
-nocall	 do not include kmer-based prediction calls into logos
-barcodelogo	Generate PWM barcode logos using colored rectangles
-heightscaledbars	Scale height of barcode logos based on nucleotide frequency
-maxheightscaledbars	Scale height of barcode logos based on maximum nucleotide frequency at position
-colorscaledbars	Scale color of barcode logos based on nucleotide frequency
-label	Include white IUPAC label to barcodelogo (cutoff = more than 50% of max)
-neg	Allow negative values in PWM (using paths)
-noname	Do not include filename in the logo
-editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)	print table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance
-longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]
-iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]
-match=[kmer1],[kmer2]	Include spacing and orientation heatmap for these kmers in kmer summary svg
-hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]
-q	print frequencies instead of counts
-n	print nucleotide counts for full sequences
-s	Print values for different gap lengths on same line
-c	Print raw counts from both files
-u	Use only unique input sequences
-p	Print p-values (Winflat program needs to be installed in path)
-i	Print incidence of all input sequences
-e	Use even background instead of background from file
-lim=[position followed by strand (one or both must be given)]	 Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)
-exc=[position][strand],[position][strand]‚..,[position][strand]	Exclude indicated positions and strands (e.g. -exc=4,3F,6R)
-m=[number]	Generate PWM using multinomial [number] distribution
-dimer	print counts for dimeric sequences only
-kmer	Print kmer count table
-dinuc	print dinucleotide data to output and svg
-ic	Output information content for all spacings
-nogaps	Count only kmers without gaps
-allgaps	Count kmers with gaps in any position (default is only middle 1 or 2 positions)
-mono	Use mononucleotide background instead of multinomial background
-both	Count both instances of palindromic hits
-bothifnotequal	Count both instances of palindromic hits if hit scores are not equal
-forwardonly	Count only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)
-reverseonly	Count only reverse instance of palindromic hits
-forwardifequal	Count only forward instance of palindromic hits if scores equal (default alternates between strands)
-reverseifequal	Count only reverse instance of palindromic hits if scores equal
-14N	Set sequence length to 14N (default is 20N)
-30N	Set sequence length to 30N
-40N	Set sequence length to 40N
-fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)
-mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits
-kl=[kmer length] kmer length used to calculate lambda (default 8)
-x	Extended output for debugging
-o=[output file name] name of svg output file



Autoseed usage:
$ ./totalautoseed -40N <Background sequence> <Signal sequence> <huddinge_distance> 
 <min_seed_length> <max_seed_length> <length_cutoff> <iupac_cutoff> <local_max_cutoff_count> <logo_number_cutoff>

Ex) In case 20N SELEX oligos.
$ ./totalautoseed -20N TTCTTC20NTTCKH3_sig.seq TTCTTC20NTTCKH4_sig.seq 1 8 10 0.35 - 50 100
TTCTTC20NTTCKH3_sig.seq and TTCTTC20NTTCKH4_sig.seq are 20 base SELEX sequence files (without headers).

Ex) In case 40N SELEX oligos.
$ ./totalautoseed -40N aaa.txt bbb.txt 1 8 10 0.35 - 50 100
aaa.txt and bbb.txt are 40 base SELEX sequence files (without headers).

Example output file:
TTCTTC20NTTCKH4_sig.seq_1_8_10_0.35_-_50_logos.svg


---
Dinucleotide heatmap generation usage:
$ ./spacek40 --f -dinuc -nocall -m=<multinomial> <Nlength_option> <Background sequence> <Signal sequence> <Seed> <cutoff>

Ex) In case 20N SELEX oligos.
$ ./spacek40 --f -dinuc -nocall -m=1 -20N TTCTTC20NTTCKH3_sig.seq TTCTTC20NTTCKH4_sig.seq NATGACGTSATN 200
TTCTTC20NTTCKH3_sig.seq and TTCTTC20NTTCKH4_sig.seq are 20 base SELEX sequence files (without headers).

Ex) In case 40N SELEX oligos.
$ ./spacek40 --f -dinuc -nocall -m=1 -40N aaa.txt bbb.txt NCCGGAAN 200
aaa.txt and bbb.txt are 40 base SELEX sequence files (without headers).

Example output file:
TTCTTC20NTTCKH4_sig.seq_logo.svg


