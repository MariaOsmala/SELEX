#!/bin/bash
#SBATCH --job-name=miss_count_kmers
#SBATCH --output=miss_count_kmers_out/count_kmers_%A_%a.out
#SBATCH --error=miss_count_kmers_err/count_kmers_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small #longrun
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00 #some need 7 days
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --array=0-22 #0-359 #3595

#35320036_[0-22] 

#sacct --format JobID%-20,State -j 31651015

#31106315_[0-10] #OK
#31106342_11-20 #OK
#31651015_[21-100] #OK
#31651026_101-200 #OK
#Time run out: 
#117,147,148,150,183,184,185
#31651038_[201-237] #OK
#Time run out: 
#223,225,226,227,237

#117,147,148,150,183,184,185,223,225,226,227,237 
#rerun: 31710601

#Still time runts out 
#31710601_150            TIMEOUT 
#31710601_227            TIMEOUT 

#rerun individually
#1500-1509 31776556 OK

#2270-2279 31776625 OK, 2271 TIMEOUT

#2238-2500: 31777882 OK, 271 TIMEOUT
#31778101_[501-750] 2501-2750 OK

#31778871_751-1000 2751-3000 #OK

#31785304_1-200 3001-3200 #OK
#31785412_201-400 3201-4300 #OK
# 31785414_401-572 3401-3572 OK

#needs longrun
# 2271,3435,3436,3438,3439,3442,3519
#rerun
# 2271: 31860440 took 5 days
# 3435,3436,3438,3439,3442,3519: 31860682 OK, took 4 days
#SBATCH --gres=nvme:5

srun count_kmers_with_Ns.sh ${SLURM_ARRAY_TASK_ID} 0


seff $SLURM_JOBID
