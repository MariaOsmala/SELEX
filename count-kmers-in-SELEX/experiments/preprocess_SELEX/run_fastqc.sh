#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastq_outs/fastqc_%A_%a.out
#SBATCH --error=fastq_errs/fastqc_%A_%a.err
#SBATCH --account=project_2013895 
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --array=301-455 # #0-455

# sacct --format JobID%-20,State -j 34048608,34048713,34048952

#34048608_0-100
#34048713_101-300
#34048952_301-344
jid="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
trap 'seff "${jid:-$SLURM_JOBID}"' EXIT

srun fastqc.sh ${SLURM_ARRAY_TASK_ID} 


# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
