cd /projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/ucsc-oligomatch
mkdir -p conda_envs/oligomatch

module load tykky

conda-containerize new --mamba --prefix conda_envs/oligomatch oligomatch.yml

export PATH="/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/ucsc-oligomatch/conda_envs/oligomatch/bin:$PATH"