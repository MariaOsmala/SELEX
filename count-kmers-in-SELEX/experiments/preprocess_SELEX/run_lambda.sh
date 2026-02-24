#!/bin/bash
#SBATCH --job-name=lambda
#SBATCH --output=outs/lambda_%A_%a.out
#SBATCH --error=errs/lambda_%A_%a.err
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10G #
#SBATCH --cpus-per-task=1
#SBATCH --array=0 # #0-357

# sacct --format JobID%-20,State -j 30097036


srun lambda.sh ${SLURM_ARRAY_TASK_ID} 

#32073986_[0] 

#SBATCH --gres=nvme:50



seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
