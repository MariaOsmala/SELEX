#!/bin/bash
#SBATCH --job-name=lambda
#SBATCH --output=lambda_outs/lambda.out
#SBATCH --error=lambda_errs/lambda.err
#SBATCH --account=project_2016851
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=1G #
#SBATCH --cpus-per-task=1


# sacct --format JobID%-20,State -j 30097036
 
##SBATCH --array=0 # #0-357

srun lambda.sh ${SLURM_ARRAY_TASK_ID} 

#32073986_[0] 

#SBATCH --gres=nvme:50



seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 3491
#svg: 2980

