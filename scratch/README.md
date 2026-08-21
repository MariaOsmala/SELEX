## Generate k-mers from seed and score them based on motif, add rank
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