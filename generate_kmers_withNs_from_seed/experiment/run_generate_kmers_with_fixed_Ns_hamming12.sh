#!/bin/bash
#SBATCH --job-name=kmers
#SBATCH --output=kmers_outs/kmers_%A_%a.out
#SBATCH --error=kmers_errs/kmers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:10
#SBATCH --array=301-364# 1-364

#sacct --format JobID%-20,Submit,Cluster,State,Elapsed,ElapsedRaw,ExitCode,User,Group,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,ReqMem,MaxRSS,TIMELIMIT,ElapsedRaw,CPUTime,TotalCPU -j 30362167,30362189,30362272,30362467,30362880 | grep FAILED



#30362167_1-10 Ok
#30362189_11-100 OK
#30362272_101-200 OK
#30362467_201-300 OK
#30362880_301-364 #OK

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
srun apptainer_wrapper exec Rscript --no-save ../code/generate_kmers_with_fixed_Ns_hamming12.R $SLURM_ARRAY_TASK_ID 0


