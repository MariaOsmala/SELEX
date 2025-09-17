#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=outs/combine.out
#SBATCH --error=errs/combine.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=360G #360
#SBATCH --cpus-per-task=1

# 28676191
# Load r-env
module load r-env/442

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2013895/tmp///" >> ~/.Renviron


# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/combine_protein_sequences_and_kmers.R $SLURM_ARRAY_TASK_ID

seff $SLURM_JOBID
