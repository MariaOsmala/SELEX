#!/bin/bash

KMERS_FILE="$1"
FASTQ_FILE="$2"

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <kmers_file> <fastq_file>"
    exit 1
fi

# Handle compressed files
read_fastq() {
    if [[ "$1" == *.gz ]]; then
        zcat "$1"
    else
        cat "$1"
    fi
}

# Create a function for parallel to use
count_in_chunk() {
    awk '
    BEGIN {
        # Load k-mers and patterns
        while ((getline kmer < "'"$KMERS_FILE"'") > 0) {
            pattern = kmer
            gsub(/N/, "[ATCG]", pattern)
            patterns[kmer] = pattern
            counts[kmer] = 0
        }
    }
    {
        # Check each sequence against all patterns
        for (kmer in patterns) {
            if ($0 ~ patterns[kmer]) {
                counts[kmer]++
            }
        }
    }
    END {
        # Output counts for this chunk
        for (kmer in patterns) {
            print kmer "\t" counts[kmer]
        }
    }'
}
export -f count_in_chunk
export KMERS_FILE

# Process in parallel and sum results
read_fastq "$FASTQ_FILE" | sed -n '2~4p' | \
parallel --pipe --block 100M -j $(nproc) count_in_chunk | \
awk '{counts[$1]+=$2} END {for(k in counts) print k"\t"counts[k]}'