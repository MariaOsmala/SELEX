#!/bin/bash
#SBATCH --job-name=run_jellyfish
#SBATCH --output=jellyfish_outs/count_8mers_%A_%a.out
#SBATCH --error=jellyfish_errs/count_8mers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=500M #
#SBATCH --cpus-per-task=4
#SBATCH --array=201-455 # #0-455

# sacct --format JobID%-20,State -j 34079075,34079103,34079262,34079465

#34083143_0-200
#34083345_201-455
srun count_jellyfish_8mers.sh ${SLURM_ARRAY_TASK_ID} 

#SBATCH --gres=nvme:50

seff $SLURM_JOBID

# dry run: list zero-byte files
# find . -type f -size 0 -print

# actually delete
# find . -type f -size 0 -delete
#kmer_counts: 2189
#lambda: 2982 
#svg: 2980
