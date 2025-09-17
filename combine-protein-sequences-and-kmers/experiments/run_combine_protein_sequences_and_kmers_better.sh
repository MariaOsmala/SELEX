#!/bin/bash
#SBATCH --job-name=better_combine
#SBATCH --output=outs/better_combine.out
#SBATCH --error=errs/better_combine.err
#SBATCH --account=project_2013895
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=5G #360
#SBATCH --cpus-per-task=1

# 28703620
# Load r-env
module load r-env/442

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2013895/tmp///" >> ~/.Renviron


# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/combine_protein_sequences_and_kmers_better.R $SLURM_ARRAY_TASK_ID

seff $SLURM_JOBID

#cd /scratch/project_2013895/SELEX/combined_kmers_protein_sequences/
#gzip -c combined_kmers_with_metadata_better_maxscore_scaled.tsv > from_sequence_to_affinity_data.tsv.gz

#Copy to allas

module load allas
allas-conf project_2013895


#a-publish -b from_sequence_to_affinity from_sequence_to_affinity_data.tsv.gz


# for file in /scratch/project_2013895/SELEX/streamed_kmers/*.tsv; do
#   if awk '{ if (NF != 3) { exit 1 } }' "$file"; then
#     : # all good
#   else
#     echo "Bad format: $file"
#   fi
# done

# These files did not have 3 columns, remove (19)
# ALX4_EOMES_CAP-SELEX_TGCACG40NTTG_AAD_NGYGYTAAYNNNNNNTNACACNN_1_3.tsv
# ALX4_TBX21_CAP-SELEX_TGGCAC40NAAG_AAD_RGGTGNTAATNNNNNNNNNCASYNN_1_3.tsv
# E2F3_EOMES_CAP-SELEX_TCTTCG40NGAT_AAA_NNGYGNNNNGGCGCSNNNNCRCNN_1_3.tsv
# ELK1_EOMES_CAP-SELEX_TGAGTC40NAAC_AAA_RGGTGNGANNNNNNNNTNNCACCGGAAGY_2_2.tsv
# ELK1_HOXA3_CAP-SELEX_TCTGTG40NAGA_AAA_NCCGGWNNNNNNNNNNNSCATTAN_1_2.tsv
# ELK1_TBX21_CAP-SELEX_TCCTAT40NAGA_AAA_RAGGTSRNNNNNNNNNNNNNNNNCGGAAGYN_2_2.tsv
# ERF_EOMES_CAP-SELEX_TGAGTC40NAAC_AAC_RRGTGTKNNNNNNNNNNNNNNCMGGANNN_1_3.tsv
# ETV2_EOMES_CAP-SELEX_TCAGCA40NCTC_AAA_NGGTGTNNNNNNNNNNNNNNNCCGGAWNNN_1_3.tsv
# ETV2_TBX21_CAP-SELEX_TGTAGC40NCCT_AAA_RGTGTKRNNNNNNNNNNNCNCMGGAARN_1_3.tsv
# ETV5_EOMES_CAP-SELEX_TACAAG40NTCC_AY_RNGTGNNNNNNNNNNNNNRCRCCGGAWSN_1_3.tsv
# ETV5_HOXA2_CAP-SELEX_TCGGCG40NACG_AY_RSCGGWAATKNNNNNNNNMATTA_2_3.tsv
# FOXJ3_TBX21_CAP-SELEX_TCCTCC40NCTC_AAE_NNGYGNNNNNNNNWAACAACACNN_1_3.tsv
# FOXO1_ETV7_CAP-SELEX_TTGGTG40NTAC_AS_RWMAACAGGNNNNNNTTCCNN_1_2.tsv
# GCM1_ETV7_CAP-SELEX_TTGGTG40NTAC_AU_NTNNNGGCGGAAGNNNTTCCNNN_2_3.tsv
# GCM1_SPDEF_CAP-SELEX_TTCTAC40NCGA_AU_RTRSKGGCGGANNNNNATCCNNN_2_3u.tsv
# GCM1_SPDEF_CAP-SELEX_TTCTAC40NCGA_AU_RTRSKGGCGGANNNNNNATCCNNN_2_3u.tsv
# HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_AGGYGYGANNNNNNNNNNNNNNNTCRTWAA_2_3.tsv
# HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_NNNACGANNNNNNTCGTNNN_1_3u.tsv
# HOXD12_TBX21_CAP-SELEX_TCCTAT40NAGA_AAB_NGGTGTNNNNNNNNNNNNNCACNTNNTWAN_1_2.tsv


