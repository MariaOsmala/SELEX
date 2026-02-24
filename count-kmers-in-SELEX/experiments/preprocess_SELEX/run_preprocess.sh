#!/bin/bash
#SBATCH --job-name=prepros_SELEX
#SBATCH --output=outs/prepros_SELEX_%A_%a.out
#SBATCH --error=errs/prepros_SELEX_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=350G #
#SBATCH --cpus-per-task=1
#SBATCH --array=252 # #0-357

# sacct --format JobID%-20,State -j 30097036

# 29906959_[0-10] #50G 30min 
# 29907097_[11-200] #50G 30min
# 29907317_201-300
#memory runs out with 50G and 100G
#30139090_218-230 #200G 12h #OK

# 30141067_231-270 200G memory run out
# 252-270 memory run out
# 30143832_252-270 300G FAILED why?
# 30167050_252


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
