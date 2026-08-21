#!/bin/bash
#SBATCH --job-name=scale
#SBATCH --output=scale_outs/scale_%A_%a.out
#SBATCH --error=scale_errs/scale_%A_%a.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:5
#SBATCH --array=101-200 #

#30220375_0-10 #OK
#30220511_11-100 #OK

#srun scale_match_scores.sh ${SLURM_ARRAY_TASK_ID}
#srun divide_match_scores_by_max.sh ${SLURM_ARRAY_TASK_ID}
#srun divide_match_scores_by_max_Hamming12.sh ${SLURM_ARRAY_TASK_ID}
srun claude_divide_scores_by_max.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID
