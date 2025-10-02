First, try to find the corresponding SELEX signal and background data for each motif and collect all data: 

correspondence-between-motifs-and-SELEX-data/code/collect_all.R

There are 3636 non Methyl-HT-SELEX motifs, for 3573 (98.3%) of these both signal 
and background are found, missing signal, background or both for 63 (1.7%). 
Missing are listed in /projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data.tsv

Try to find the protein sequences for each motif

/projappl/project_2013895/SELEX/TF-amino-acid-sequences/

code/ #scripts
Data/ #Supplementary tables

#try to extract protein sequences from supplementary info
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/code/add_protein_sequence.R


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