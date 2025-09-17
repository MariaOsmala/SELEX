
#3360
mapfile -t scaled_kmers_for_included_motifs < /projappl/project_2013895/SELEX/combine-protein-sequences-and-kmers/experiments/included_motifs.txt


cd /scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore #3733

output=/scratch/project_2013895/SELEX/streamed_kmers_scaled_by_maxscore_included
mkdir $output

for f in "${scaled_kmers_for_included_motifs[@]}"; do
  cp $f".tsv" "$output/"
done

cd $output

for f in "${scaled_kmers_for_included_motifs[@]}"; do
  gzip $f".tsv"
done

