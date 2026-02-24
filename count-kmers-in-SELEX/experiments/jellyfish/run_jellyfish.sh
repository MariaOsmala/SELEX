#!/bin/bash
#SBATCH --job-name=jellyfish
#SBATCH --output=outs/jellyfish_%A_%a.out
#SBATCH --error=errs/jellyfish_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --array=201-357 # #0-357

# sacct --format JobID%-20,State -j 30097036

#filtered data:
#30354552_0-10 #OK
#30354568_11-200 #OK
#30354758_201-357

#non-filtered data:
# 30144920_0 OK 30min 5G
# 30144938_1-100 OK 30min 5G 
# 30146719_101-200 12h OK
# 30166582_201-300 #OK
#  30167074_301-357 #OK
srun jellyfish.sh ${SLURM_ARRAY_TASK_ID} 

#SBATCH --gres=nvme:50

seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
