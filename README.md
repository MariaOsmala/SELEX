#README

## Find the corresponding SELEX signal and background data for each motif and collect all data:

```{sh}
correspondence-between-motifs-and-SELEX-data/code/collect_all.R
```

There are 3635 non Methyl-HT-SELEX motifs, for 3596 (98.9%) of these both signal
and background are found, missing signal, background or both for 39 (1.1%).
Missing are listed in `/projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data_19032026.tsv`

Motifs with SELEX data are in `motifs_with_SELEX_signal_and_background_19032026.tsv` 3597 rows

Metadata for all 3933 motifs in `metadata_final_with_SELEX_data_19032026.tsv` Dated 2.4.2026

## Try to find the protein sequences for each motif
```{sh}
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/

code/ #scripts
Data/ #Supplementary tables
```

Try to extract protein sequences from supplementary info
```{sh}
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/code/add_protein_sequence.R
```
Create "TF_info_Ensembl115.txt"

Lists all TFs for which there is a motif:
`TF-info/code/TF_interaction_summary.R` `creates RData/TFs_with_motifs.RDS`
` TF-info/code/TF_alternative_names.R` `creates TF_info_Ensembl115.txt`

Check protein sequences
```{sh}
TF-amino-acid-sequences/code/check_protein_sequences_all_motifs.R
```
Protein sequence is missing or wrong for 28 motifs

## Combine motifs, SELEX data and protein sequences
/projappl/project_2013895/SELEX/TF-amino-acid-sequences/code/combine_motifs_SELEXdata_proteinSequences.R

Output in "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_and_proteins_01042026.tsv"
contains also the links to downloadable selex files

## Load SELEX data to allas:

/projappl/project_2013895/SELEX/data-to-allas/data_to_allas.sh

##Generate k-mers 

Generate canonical realisations and their Hamming1 and Hamming2 neighborhood. 
Keep Ns as N and generate substitutions only at positions that are not N.

```{sh}
/SELEX/generate_kmers_withNs_from_seed/experiment/run_generate_kmers_with_fixed_Ns_hamming12.sh
/SELEX/generate_kmers_withNs_from_seed/code/generate_kmers_with_fixed_Ns_hamming12.R
```
Generated canonical k-mers are in `/scratch/project_2013895/SELEX/kmers_from_fixed_N_seeds/`
Separately for realisations, hamming1, and  hamming2.

Combine k-mers, Hamming1 or Hamming1 + Hamming2 neighborhood

```{sh}
SELEX/generate_kmers_withNs_from_seed/experiment/run_combine_kmers.sh
```

Results are as ".txt" and "_hamming1.txt" in
`/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/`

## Count the occurrences of generated k-mers in SELEX data

In signal and background.

Counts reads containing the k-mer, not k-mer occurences.
Count is 1 even if a k-mer occurs in multiple positions of the same read

```{sh}
/SELEX/count-kmers-in-SELEX/experiments/count_kmers_with_Ns/run_count_kmers_with_Ns.sh
/SELEX/count-kmers-in-SELEX/experiments/count_kmers_with_Ns/count_kmers_with_Ns.sh
```
Output in `/scratch/project_2013895/SELEX/kmer_counts_combined_fixed_N_seeds/`

Counts all occurrences of k-mers. 
TODO

#Score k-mers, subtract background from signal, consider the non-specific carry over
/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/experiment/run_score_kmers.sh
/projappl/project_2013895/SELEX/generate_kmers_withNs_from_seed/code/match_motifs_to_kmers.R

results are in /scratch/project_2013895/SELEX/scored_kmers_fixed_N_seeds_June2026
tables with: Kmer    Corrected_count Corrected_count_scaled  Corrected_count_rank    motif_match_score       motif_match_score_scaled        motif_match_score_rank

## SELEX data preprocessing

There are 5442 unique SELEX signal or background files. Nine signal .fastq files are corrupted.

Check `SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/code/unique_SELEX_data_files.R`

Run fastqc:
```{sh}
/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/run_fastqc.sh
```
The output from fastqc analysis are in `/projappl/project_2013895/SELEX/quality_control/fastqc`

Remove reads with N
`/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/run_seqkit.sh`

## Lambda computation

There are 2681 unique combinations of SELEX signal and data. 2672 if corrupted .fastq files are removed. 
`SELEX/count-kmers-in-SELEX/experiments/lambda_computation/code/unique_combinations_of_signal_and_background.R`

Compute lambda using spacek. This converts .fastq -> .fasta -> .seq. Runs happily with corrupted signal files, that are
probably corrected in the conversion. 
`/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/run_lambda.sh`

My own implementation of lambda computation.
`/projappl/project_2013895/SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments/run_compute_lambda_R.sh`

Spacek and my own implementation give the same result, check the following code. Draws also lambda distribution.
`SELEX/count-kmers-in-SELEX/experiments/preprocess_SELEX/code/parse_lambda.R`

Experimented also with jellyfish. This seems to result to the same lambda values.
`SELEX/count-kmers-in-SELEX/experiments/lambda_computation/experiments/run_jellyfish_count8mers.sh`
`SELEX/count-kmers-in-SELEX/experiments/lambda_computation/code/compute_lambda_jellyfish.R`

# For which motifs the score computation was successfull
generate_kmers_withNs_from_seed/code/lambda_distribution.R
generates: "/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026_lambda_info.tsv",
# Protein sequences and stoichiometry
TF-amino-acid-sequences/code/combine_motifs_SELEXdata_proteinSequences.R

# Combine protein sequences, stoichiometries and k-mer counts
/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/code/combine_protein_sequences_and_kmers_better.R
/projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/run_combine_protein_sequences_and_kmers_better.sh

#Data will be located in 
/scratch/project_2013895/SELEX/sequence-to-affinity-data-01072026

all.tsv 7.8GB

ls separate_with_protein_sequences/ | wc -l 3277
ls separate_without_protein_sequences/ | wc -l 3277

#make tar.gz
cd /scratch/project_2013895/SELEX
# tar -czf sequence-to-affinity-data-01072026.tar.gz sequence-to-affinity-data-01072026/

#make zip

zip -r -X -9 sequence-to-affinity-data-01072026.zip sequence-to-affinity-data-01072026/ 