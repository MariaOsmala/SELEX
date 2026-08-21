#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=combine_outs/combine_%A_%a.out
#SBATCH --error=combine_errs/combine_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:5
#SBATCH --array=11-354 #354 

# 30200892_0-10
# 30200915_11-354
srun claude_combine_kmers.sh ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
