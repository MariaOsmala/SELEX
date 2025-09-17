#!/bin/bash
#SBATCH --job-name=scale
#SBATCH --output=scale_outs/scale_%A_%a.out
#SBATCH --error=scale_errs/scale_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:5
#SBATCH --array=11-96 #96 #346 -377 # #0-377 #3774  



#29402165_0-10 OK
#29402299_11-200 OK
# 29405872_201-346 OK

#OLD
#28702155_0-10 OK
#28702361_11-100 Ok
#28702552_101 #OK
# 28702549_102-377




#srun scale_match_scores.sh ${SLURM_ARRAY_TASK_ID}
#srun divide_match_scores_by_max.sh ${SLURM_ARRAY_TASK_ID}
srun divide_match_scores_by_max_Hamming12.sh ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
