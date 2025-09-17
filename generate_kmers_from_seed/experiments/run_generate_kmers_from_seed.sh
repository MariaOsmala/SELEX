#!/bin/bash
#SBATCH --job-name=kmers
#SBATCH --output=outs/kmers_%A_%a.out
#SBATCH --error=errs/kmers_%A_%a.err
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=50G #
#SBATCH --cpus-per-task=1
#SBATCH --array=601-933 #3933


# 28590553
#sacct --format JobID%-20,State -j 28590553 | grep FAILED

#28593088 1-300 #OK
#Max memory: 1G 
#Max elapsed time (seconds): 1806=30 min

#28593089 301-600 #OK 1G 47 min
#28593090 601-1000 #1G 22H
# Failes due long-gapped seed
# 821,822,834,838,847,,857,861,876,895,898,906,920,951,991,992,999

# 28607110 1001-1300

# 18,81,127,136,138,144,145,150,153,174,176,194,205,206,209,212,216,220,223,224,238,
# 239,244,246,247,248,250,254,270,276,281,282,285,286,288,289,292,293,294,295,296,297 

# 28607134 1301-1600, 1990
# 305,306,310,311,314,315,316,317,325,326,345,365,374,375,376,377,424

#28602804 1601-1900 #OK
#28602968 1901-2000 OK
# 

#28613719 2001-2300
# 254,278,279 

#28613821 2301-2600
# 305,306,376,393,557 

#28613851 2601-2900
#607,615,616,646,648,650,660,682,684,694,705,718,719,733,734,749,772,831,832 



#28613854 2901-3000 #OK

#28620125 3001-3300
#279, 283

#28620142 3301-3600
#425,426,397

#28620756 3601-3933 #one from this still running
#731,906


# sacct --format=JobID,MaxRSS,Elapsed -j 28593090 --parsable2 | \
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
