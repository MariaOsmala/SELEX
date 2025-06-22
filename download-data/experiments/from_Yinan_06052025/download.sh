#!/bin/bash
#SBATCH --job-name=download
#SBATCH --output=outs/download.out
#SBATCH --error=errs/download.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2G

cd /scratch/project_2013895/SELEX/data/
mkdir Yinan_M0_M33_M67_M100_06052025
cd Yinan_M0_M33_M67_M100_06052025

# Tried in several ways but did not work. How to download data from funet 
#file server directly to csc (if requires loggin in to the system)

wget "https://filesender.funet.fi/?s=download&token=025a21f7-5f29-451f-9f2a-9c181f64f694" -O data_from_FS
wget https://filesender.funet.fi/download.php?files_ids=1030347

https://filesender.funet.fi/?s=download&token=025a21f7-5f29-451f-9f2a-9c181f64f694


#56G mCSELEX_fastq_trimmed.zip
#59 393 571 793 mCSELEX_fastq_trimmed.zip