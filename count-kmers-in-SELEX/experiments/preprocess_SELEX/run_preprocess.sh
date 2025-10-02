#!/bin/bash
#SBATCH --job-name=prepros_SELEX
#SBATCH --output=outs/prepros_SELEX_%A_%a.out
#SBATCH --error=errs/prepros_SELEX_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=50G #
#SBATCH --cpus-per-task=1
#SBATCH --array=301-357 # #0-357

# 29906959_[0-10]
# 29907097_[11-200]
# 29907317_201-300
# 29910498_301-357
srun preprocess.sh ${SLURM_ARRAY_TASK_ID} 

#SBATCH --gres=nvme:50

seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
