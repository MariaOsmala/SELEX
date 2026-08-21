#!/bin/bash
#SBATCH --job-name=better_kmers
#SBATCH --output=outs/better_kmers_%A_%a.out
#SBATCH --error=errs/better_kmers_%A_%a.err
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:10
#SBATCH --array=901-933 #1-3933

#sacct --format JobID%-20,State -j 28590553 | grep FAILED

#28635937 1-10 OK
#28640369 11-300 OK
#28640380 301-600 #OK
#28640384 601-900 #OK
#28640590 901-1000 #OK

#991 Memory run out 
#999

#Time run out
#906,920,941,950,951,992,931 


#28650165 1001-1300 
#TIMEOUT
#18,81,127,136,144,145,150,153,157,174,176,177,183,193,194,205,206,209,212,216,220,223,224,237,238,239,240,243,244,246,247,248,250,254,255,270,271,276
#282,285,286,287,288,289,291,292,294,295,296,297

#Memory runs out 
#293,281,138

#FAILED
#36,46,147,152,163,165,166,170,175,178,182,185,190,192,203,207,219,229,241,256,264,290,299 



#28650167 1301-1600 
#TIME RUNS OUT
#305,307,310,311,313,314,315,316,317,325,326,329,344,345,368,369,374,375,377

#MEMORY RUNS OUT
#306,365,376 

#301,312,324,330,332,334,342,343,370,424,527


#28650169 1601-1900 #OK
#28650213 1901-2000 #OK
#990            TIMEOUT

#28657013 2001-2300 #
#TIMEOUT
#221,230,254,278,279,296,

#FAILED
#76,153,286,290

#28657014 2301-2600 #
#305,306,311,376,393,424,444,557,579,308,401,502,556,559,561,589

#28657021 2601-2900 #
#607,615,616,618,635,643,646,648,650,660,665,682,684,694,695,699,705,716,718 
#719,726,733,734,749,762,772,831,832 

#640,642,675,683,717,725,736,742 


#28657085 2900-3000 #OK

#28662575 3001-3300 #
# 279,283,296,76,77,78,141

#28662577 3301-3600 #some still running
#28662580 3600-3900 #
#606,679,682,731,657,742,746,884 
#28662584 3901-3933 #
#905,906,918


# sacct --format=JobID,MaxRSS,Elapsed -j 28695827 --parsable2 | \
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
srun apptainer_wrapper exec Rscript --no-save ../code/generate_kmers_from_seed_better.R $SLURM_ARRAY_TASK_ID




seff $SLURM_JOBID
