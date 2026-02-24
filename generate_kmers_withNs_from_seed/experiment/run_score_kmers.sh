#!/bin/bash
#SBATCH --job-name=score_kmers
#SBATCH --output=score_kmers_outs/score_kmers_%A_%a.out
#SBATCH --error=score_kmers_errs/score_kmers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --array=1-358 # 1-358

#32031018_[1-358] 

#sacct --format JobID%-20,State -j 30043694 | grep FAILED

#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node> #SBATCH --gres=nvme:50
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

# Load r-env
module load r-env/442

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2013895/tmp///" >> ~/.Renviron


jid="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
trap 'seff "${jid:-$SLURM_JOBID}"' EXIT

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/match_motifs_to_kmers.R $SLURM_ARRAY_TASK_ID 0


