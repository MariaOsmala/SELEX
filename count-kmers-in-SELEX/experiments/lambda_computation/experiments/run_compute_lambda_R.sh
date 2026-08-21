#!/bin/bash
#SBATCH --job-name=run_lambda
#SBATCH --output=lambda_outs/%A_%a.out
#SBATCH --error=lambda_errs/%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --array=0-267 # #0-267

# sacct --format JobID%-20,State -j 34126530, 34126551, 34126879

# 34130334

#run first: /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/code/unique_combinations_of_signal_and_background.R

srun compute_lambda_R.sh ${SLURM_ARRAY_TASK_ID} 

#SBATCH --gres=nvme:50

seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
