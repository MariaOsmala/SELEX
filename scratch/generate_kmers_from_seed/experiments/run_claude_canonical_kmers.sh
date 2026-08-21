#!/bin/bash
#SBATCH --job-name=claude_kmers
#SBATCH --output=kmers_outs/claude_kmers_%A_%a.out
#SBATCH --error=kmers_errs/claude_kmers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=100G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:30
#SBATCH --array=238,395-396,398,412,454,467,469,472-478,482,486,489-490,492,512-514,517-518,522-523,526,531-532,535-536,539-545,547-551,554-557,559,561-563,567-585,587-590,592-594,597-602,604,606-607,609-635  # # #91-635 

#sacct --format JobID%-20,State -j 30043694 | grep FAILED

# 30043694 1-10
# 30043716 11-100
# 30043813 101-200
# 30043923 201-250
# 30044187 251-300
# 30044415 301-364 successfull until 309

# 30046865_91-200 3091-3200 OK
# 30047000_201-400 OK

# 30047621_401-635

#Failed 412

#Memory runs out
# rerun 238 395 
# 30056742_238,395 50G

#238,395-396,398,412,454,467,469,472-478,482,486,489-490,492,512-514,517-518,522-523,526,531-532,535-536,539-545,547-551,554-557,559,561-563,567-585,587-590,592-594,597-602,604,606-607,609-635
#30076006 100G

#TIMEOUT


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
srun apptainer_wrapper exec Rscript --no-save ../code/claude_canonical_kmers.R $SLURM_ARRAY_TASK_ID 3000


