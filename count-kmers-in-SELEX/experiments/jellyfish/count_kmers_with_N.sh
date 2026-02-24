#!/bin/bash

KMERS_FILE="$1"        # File with k-mers containing Ns
COUNTS_FILE="$2"       # File with k-mer counts from jellyfish

# Function to get reverse complement
reverse_complement() {
    echo "$1" | tr 'ATCGN' 'TAGCN' | rev
}

echo -e "K-mer\tCount"

while read kmer_pattern; do
    [ -z "$kmer_pattern" ] && continue
    
    # Get reverse complement of the pattern
    rc_pattern=$(reverse_complement "$kmer_pattern")
    
    # Convert N to . for regex
    forward_regex=$(echo "$kmer_pattern" | sed 's/N/./g')
    reverse_regex=$(echo "$rc_pattern" | sed 's/N/./g')
    
    # Match both forward and reverse patterns
    # (avoid double counting if pattern is self-complementary)
    if [ "$forward_regex" == "$reverse_regex" ]; then
        # Self-complementary pattern
        count=$(awk -v pat="^${forward_regex}$" \
                '$1 ~ pat {sum += $2} END {print sum+0}' "$COUNTS_FILE")
    else
        # Match both orientations
        count=$(awk -v fpat="^${forward_regex}$" -v rpat="^${reverse_regex}$" \
                '$1 ~ fpat || $1 ~ rpat {sum += $2} END {print sum+0}' "$COUNTS_FILE")
    fi
    
    echo -e "$kmer_pattern\t$count"
done < "$KMERS_FILE"