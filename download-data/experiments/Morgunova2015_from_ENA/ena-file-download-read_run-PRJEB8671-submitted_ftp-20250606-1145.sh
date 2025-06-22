#!/bin/bash
#SBATCH --job-name=sub_download
#SBATCH --output=outs/sub_download.out
#SBATCH --error=errs/sub_download.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G


# Project: PRJEB8671
# The mammalian cell cycle is controlled by the E2F family of transcription factors. 
# Typical E2Fs bind to DNA as heterodimers with the related DP proteins, 
# whereas the atypical E2Fs, E2F7 and E2F8, contain two DNA-binding domains 
# (DBDs) and act as repressors. To understand the mechanism of repression we have
# resolved the structure of E2F8 in complex with DNA at atomic resolution. 
# We find that the first and second DBDs of E2F8 resemble the DBDs of typical E2F
# and DP proteins, respectively. Using molecular dynamics simulations, 
# biochemical affinity measurements and chromatin immunoprecipitation, 
# we further show that both atypical and typical E2Fs bind to similar DNA-sequences 
# in vitro and in vivo. Our results represent the first crystal structure of an 
# E2F protein with two DBDs, and reveal the mechanism by which atypical E2Fs 
# can repress canonical E2F target genes and exert their negative influence 
# on cell cycle progression.

dir=/scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp

mkdir -p $dir

cd $dir


wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074952/E2F2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074951/TFDP1.fastq.gz
