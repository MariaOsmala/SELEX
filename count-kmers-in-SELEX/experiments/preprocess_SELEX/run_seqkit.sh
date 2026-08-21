#!/bin/bash
#SBATCH --job-name=seqkit
#SBATCH --output=seqkit_outs/seqkit_%A_%a.out
#SBATCH --error=seqkit_errs/seqkit_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --array=301-455 # #0-455

# sacct --format JobID%-20,State -j 34047706,34047750,34047941,34048003

#34047706_0-10
#34047750_11-200
#34047941_201-300
#34048003_301-455
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
