#!/bin/bash
array=$1 #0-455

module load fastqc/0.12.1

file="unique_SELEX.txt" #4552
#output=/scratch/project_2013895/SELEX/fastqc
output=/scratch/project_2013895/SELEX/fastqc_reads_with_N_removed

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

#The following fastq files are corrupted
#/scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/HOXD3_CREM_3_YLIII_TGTAAG40NACT.fastq.gz
#Failed to process file HOXD3_CREM_3_YLIII_TGTAAG40NACT.fastq.gz
# uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Midline 'GACTTTTAGGCAGCCCGACCTGCAAATCTGCCACGACAAT' didn't start with '+' at 811523
#         at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:179)
#         at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:129)
#         at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
#         at java.lang.Thread.run(Thread.java:748)

#zcat /scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/HOXD3_CREM_3_YLIII_TGTAAG40NACT.fastq.gz | sed -n '811521,811525p'

# /scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/HOXD4_CREM_3_YLIII_TTCTGG40NTGC.fastq.gz
# application/gzip
# Failed to process file HOXD4_CREM_3_YLIII_TTCTGG40NTGC.fastq.gz
# uk.ac.babraham.FastQC.Sequence.SequenceFormatException: ID line didn't start with '@' at line 752453
#         at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:163)
#         at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:129)
#         at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
#         at java.lang.Thread.run(Thread.java:748)

#zcat /scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp/HOXD4_CREM_3_YLIII_TTCTGG40NTGC.fastq.gz | sed -n '752451,752455p'


for ((i=start_ind; i<=end_ind; i++)); do
   echo $i
   echo ${files[$i]}
   file=${files[$i]}
   base="${file##*/}"    
   #replace the file name
   #fastqc --extract -o $output ${files[$i]} --quiet #orig files
   fastqc --extract -o $output /scratch/project_2013895/SELEX/data_exclude_reads_with_N/$base --quiet #reads with N removed
   
done



# "ELF1_CREM_3_YLIII_TTCATT40NTAT"  1747
# "ELF4_CREM_3_YLIII_TCTAGA40NCGA"  1768
# "HOXA4_CREM_3_YLIII_TAGACT40NGTG" 1987
# "HOXD3_CREM_3_YLIII_TGTAAG40NACT" Corrupted fastq
# "HOXD4_CREM_3_YLIII_TTCTGG40NTGC" Corrupted fastq
# "PAX1_VSX2_3_YLIIII_TATAGC40NGAC" 2263 
# "PAX2_VSX2_3_YLIIII_TTGCTG40NCTA" 2283 
# "POU4F3_CREM_3_YLIII_TAAGTT40NGAT" 2346
# "RFX5_CREM_3_YLIII_TCCCTT40NCAT"   2367


#FastQC - A high throughput sequence QC analysis tool

# SYNOPSIS
# 
#         fastqc seqfile1 seqfile2 .. seqfileN
#         fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
#            [-c contaminant file] seqfile1 .. seqfileN

# DESCRIPTION
# 
#     FastQC reads a set of sequence files and produces from each one a quality
#     control report consisting of a number of different modules, each one of 
#     which will help to identify a different potential type of problem in your
#     data.
# 
#    If no files to process are specified on the command line then the program
#     will start as an interactive graphical application.  If files are provided
#     on the command line then the program will run with no user interaction
#     required.  In this mode it is suitable for inclusion into a standardised
#     analysis pipeline.
#     
#     The options for the program as as follows:
#     
#     -h --help       Print this help file and exit
#     
#     -v --version    Print the version of the program and exit
#     
#     -o --outdir     Create all output files in the specified output directory.
#                     Please note that this directory must exist as the program
#                     will not create it.  If this option is not set then the 
#                     output file for each sequence file is created in the same
#                     directory as the sequence file which was processed.
#                     
#     --casava        Files come from raw casava output. Files in the same sample
#                     group (differing only by the group number) will be analysed
#                     as a set rather than individually. Sequences with the filter
#                     flag set in the header will be excluded from the analysis.
#                     Files must have the same names given to them by casava
#                     (including being gzipped and ending with .gz) otherwise they
#                     won't be grouped together correctly.
#                     
#     --nano          Files come from nanopore sequences and are in fast5 format. In
#                     this mode you can pass in directories to process and the program
#                     will take in all fast5 files within those directories and produce
#                     a single output file from the sequences found in all files.                    
#                     
#     --nofilter      If running with --casava then don't remove read flagged by
#                     casava as poor quality when performing the QC analysis.
#                    
#     --extract       If set then the zipped output file will be uncompressed in
#                     the same directory after it has been created. If --delete is 
#                     also specified then the zip file will be removed after the 
#                     contents are unzipped. 
#     -j --java       Provides the full path to the java binary you want to use to
#                     launch fastqc. If not supplied then java is assumed to be in
#                     your path.
#                    
#     --noextract     Do not uncompress the output file after creating it.  You
#                     should set this option if you do not wish to uncompress
#                     the output when running in non-interactive mode.
#                     
#     --nogroup       Disable grouping of bases for reads >50bp. All reports will
#                     show data for every base in the read.  WARNING: Using this
#                     option will cause fastqc to crash and burn if you use it on
#                     really long reads, and your plots may end up a ridiculous size.
#                     You have been warned!
#                     
#     --min_length    Sets an artificial lower limit on the length of the sequence
#                     to be shown in the report.  As long as you set this to a value
#                     greater or equal to your longest read length then this will be
#                     the sequence length used to create your read groups.  This can
#                     be useful for making directly comaparable statistics from 
#                     datasets with somewhat variable read lengths.
# 
#     --dup_length    Sets a length to which the sequences will be truncated when 
#                     defining them to be duplicates, affecting the duplication and
#                     overrepresented sequences plot.  This can be useful if you have
#                     long reads with higher levels of miscalls, or contamination with
#                     adapter dimers containing UMI sequences.
# 
#                     
#     -f --format     Bypasses the normal sequence file format detection and
#                     forces the program to use the specified format.  Valid
#                     formats are bam,sam,bam_mapped,sam_mapped and fastq
#                     
# 
#     --memory        Sets the base amount of memory, in Megabytes, used to process 
#                     each file.  Defaults to 512MB.  You may need to increase this if
#                     you have a file with very long sequences in it.
# 
#     --svg           Save the graphs in the report in SVG format.
# 
#     -t --threads    Specifies the number of files which can be processed
#                     simultaneously.  Each thread will be allocated 250MB of
#                     memory so you shouldn't run more threads than your
#                     available memory will cope with, and not more than
#                     6 threads on a 32 bit machine
#                   
#     -c              Specifies a non-default file which contains the list of
#     --contaminants  contaminants to screen overrepresented sequences against.
#                     The file must contain sets of named contaminants in the
#                     form name[tab]sequence.  Lines prefixed with a hash will
#                     be ignored.
# 
#     -a              Specifies a non-default file which contains the list of
#     --adapters      adapter sequences which will be explicity searched against
#                     the library. The file must contain sets of named adapters
#                     in the form name[tab]sequence.  Lines prefixed with a hash
#                     will be ignored.
#                     
#     -l              Specifies a non-default file which contains a set of criteria
#     --limits        which will be used to determine the warn/error limits for the
#                     various modules.  This file can also be used to selectively 
#                     remove some modules from the output all together.  The format
#                     needs to mirror the default limits.txt file found in the
#                     Configuration folder.
#                     
#    -k --kmers       Specifies the length of Kmer to look for in the Kmer content
#                     module. Specified Kmer length must be between 2 and 10. Default
#                     length is 7 if not specified.
#                     
#    -q --quiet       Suppress all progress messages on stdout and only report errors.
#    
#    -d --dir         Selects a directory to be used for temporary files written when
#                     generating report images. Defaults to system temp directory if
#                     not specified.
