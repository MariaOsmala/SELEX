#!/bin/bash
#SBATCH --job-name=kmers
#SBATCH --output=outs/kmers_%A_%a.out
#SBATCH --error=errs/kmers_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=50G #
#SBATCH --cpus-per-task=1
#SBATCH --array=0-393 # #0-393  


# 28590553


# sacct --format=JobID,MaxRSS,Elapsed -j 28553183 --parsable2 | \
# awk -F'|' '
#   function to_sec(t) {
#     split(t, a, ":")
#     return a[1]*3600 + a[2]*60 + a[3]
#   }
# 
#   $2 ~ /^[0-9]+K$/ {
#     mem = substr($2, 1, length($2)-1)
#     if (mem > max_mem) max_mem = mem
#   }
# 
#   $3 ~ /^[0-9]{2}:[0-9]{2}:[0-9]{2}$/ {
#     sec = to_sec($3)
#     if (sec > max_time) max_time = sec
#   }
# 
#   END {
#     printf "Max memory (KB): %d\n", max_mem
#     print "Max elapsed time (seconds):", max_time
#   }
#   '




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


# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/generate_kmers_from_seed.R $SLURM_ARRAY_TASK_ID




seff $SLURM_JOBID
