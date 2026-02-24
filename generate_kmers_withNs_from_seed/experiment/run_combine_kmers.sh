#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=combine_outs/combine_%A_%a.out
#SBATCH --error=combine_errs/combine_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:5
#SBATCH --array=0-364 #0-364 


#sacct --format JobID%-20,State -j
#30743797_0-364 

srun combine_kmers.sh ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
