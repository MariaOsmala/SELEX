#!/bin/bash
#SBATCH --job-name=score_claude_kmers
#SBATCH --output=score_kmers_outs/score_claude_kmers_%A_%a.out
#SBATCH --error=score_kmers_errs/score_claude_kmers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00 #12hours
#SBATCH --mem-per-cpu=50G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:30
#SBATCH --array=601-700  #1-364 #

#sacct --format JobID%-20,State -j 30043694 | grep FAILED

#30356621_1-10 #100G OK
#30356764_11-100 #100G OK
#30356773_101-200 #100G 189 memory run out
#rerun 30361300_189 OK
#30356778_201-230 #100G OK

#30361324_301-400 #2301-2400 150G OK
#30362121_401-500 #OK
#30363046_501-600 
# 30363046_580 and onwards failed because the SELEX data is empty
#30371892_580 OK
#30371933_581-600
#30372317_600-700 #50

#Old
#30206024_1-100 50G
# 100 failed
#30207129_101-200 100G
#FAILED 143,157,161

#rerun 30218540_100,143,157,161 #OK
# 157 still failed, try [1] 1561-1570 30220514 OK

#30207131_201-230 #100 OK
#30218482_231-270 100 FAILS start from 2301 ->
#30220835_301-400 #2301-2400

#memory run out, rerun 30221983 OK
#30220835_318.0       CANCELLED+ 
#30220835_319.0       CANCELLED+ 

#30221998_401-600 2401-2600

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
srun apptainer_wrapper exec Rscript --no-save ../code/score_claude_canonical_kmers.R $SLURM_ARRAY_TASK_ID 2000


