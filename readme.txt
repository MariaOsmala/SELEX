First, try to find the corresponding SELEX signal and background data for each motif and collect all data:

correspondence-between-motifs-and-SELEX-data/code/collect_all.R

There are 3635 non Methyl-HT-SELEX motifs, for 3596 (98.9%) of these both signal
and background are found, missing signal, background or both for 39 (1.1%).
Missing are listed in /projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data_19032026.tsv

Motifs with SELEX data are in motifs_with_SELEX_signal_and_background_19042026.tsv

Metadata for all 3933 motifs in metadata_final_with_SELEX_data_19032026.tsv

Try to find the protein sequences for each motif

/projappl/project_2013895/SELEX/TF-amino-acid-sequences/

code/ #scripts
Data/ #Supplementary tables

#try to extract protein sequences from supplementary info
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/code/add_protein_sequence.R

# create "TF_info_Ensembl115.txt"

# TF-info/code/TF_interaction_summary.R creates RData/TFs_with_motifs.RDS
# lists all TFs for which there is a motif

# TF-info/code/TF_alternative_names.R creates TF_info_Ensembl115.txt

# check protein sequences
TF-amino-acid-sequences/code/check_protein_sequences_all_motifs.R

#combine motifs, SELEX data and protein sequences
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/code/combine_motifs_SELEXdata_proteinSequences.R

output in "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_25082025.tsv"
contains also the links to downloadable selex files:

#load SELEX data to allas:

/projappl/project_2013895/SELEX/data-to-allas/data_to_allas.sh

Generate k-mers from seed and score them based on motif, add rank
/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/run_generate_kmers_from_seed_better.sh

Only for non-Methyl-SELEX motifs:

/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/run_generate_kmers_from_seed_exclude_Methyl-HT-SELEX.sh

Results are in
/scratch/project_2013895/SELEX/streamed_kmers/

Divide the match scores by max:

/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/run_scale_match_scores.sh
/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/divide_match_scores_by_max.sh

results are in /scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore/

Try to count k-mers and score them based on SELEX data
/projappl/project_2013895/SELEX/count-kmers-in-SELEX


# Generated k-mer realisations of motif seeds and their Hamming distance 1 neighbours (base substitutions only at positions that are not N )
/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/run_claude_canonical_kmers.sh
/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/claude_canonical_kmers.R

#results in /scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/

#combine k-mers
/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/run_claude_combine_kmers.sh
/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/claude_combine_kmers.sh

#Generate k-mer realisations only at positions that are not Ns

#Hamming1 and Hamming2 neighbourhood
/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/experiment/run_generate_kmers_with_fixed_Ns_hamming12.sh

#Combine k-mers, Hamming1 or Hamming1 + Hamming2 neighbourdhood
#Check that Hamming2 does not contain Hamming1?
/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/experiment/run_combine_kmers.sh

#k-mer counts in SELEX reads, reads with Ns removed

#remove reads
/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/run_seqkit.sh

#compute lambda
/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/run_lambda.sh

#Combine the data into a table:

