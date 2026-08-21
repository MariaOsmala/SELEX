#!/bin/bash
#SBATCH --job-name=seqkit
#SBATCH --output=outs/seqkit_%A_%a.out
#SBATCH --error=errs/seqkit_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --array=201-357 # #0-357

# sacct --format JobID%-20,State -j 30097036


jid="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
trap 'seff "${jid:-$SLURM_JOBID}"' EXIT

srun seqkit.sh ${SLURM_ARRAY_TASK_ID} 

#SBATCH --gres=nvme:50


# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
