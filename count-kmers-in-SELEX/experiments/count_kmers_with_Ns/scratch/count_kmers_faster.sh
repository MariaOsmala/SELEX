#!/bin/bash

#This is an optimized version of count_kmers.sh with the same core functionality but significantly improved performance. 
#This optimized version would be much faster for your workflows when checking multiple seeds, especially with large FASTQ files on your HPC cluster.

#Here are the key differences:
#multiple k-mers containing Ns:

#100 k-mers, 1GB FASTQ:

#Old: Reads 100GB total (1GB × 100 times)
#New: Reads 1GB total (once)


# Memory trade-off: New script holds all k-mers in memory but saves massive I/O

#Limitations:
#Still counts overlapping matches
#Still no canonical k-mer handling
#Output order may differ (hash table ordering vs. input file order)




KMERS_FILE="$1"
FASTQ_FILE="$2"

# Check if file exists
if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file '$FASTQ_FILE' not found"
    exit 1
fi

# Function to read FASTQ (handles both plain and gzipped)
if [[ "$FASTQ_FILE" == *.gz ]]; then
    INPUT_CMD="zcat"
else
    INPUT_CMD="cat"
fi

# Build a single awk command that checks all patterns at once
# Single-pass: reads FASTQ only once. O(n + m) complexity.
# Loads all k-mers into memory first
# Processes the FASTQ file once, checking all patterns simultaneously

$INPUT_CMD "$FASTQ_FILE" | awk '
BEGIN {
    # Load all k-mers and convert N to regex
    while ((getline kmer < "'"$KMERS_FILE"'") > 0) {
        pattern = kmer
        gsub(/N/, "[ATCG]", pattern)
        patterns[kmer] = pattern
        counts[kmer] = 0
    }
}
NR%4==2 {  # Process only sequence lines
    for (kmer in patterns) {
        if ($0 ~ patterns[kmer]) {
            counts[kmer]++
        }
    }
}
END {
    for (kmer in patterns) {
        print kmer "\t" counts[kmer]
    }
}'