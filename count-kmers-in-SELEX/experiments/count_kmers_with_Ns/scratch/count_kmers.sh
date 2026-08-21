#!/bin/bash

# The script searches for and counts specific k-mer patterns (DNA sequences) within FASTQ files, 
# treating 'N' characters as wildcards that can match any nucleotide (A, T, C, or G).

# Limitations to Consider:
# Given your workflow, this script has some important limitations:
# It counts overlapping occurrences within sequences
# It doesn't handle canonical k-mers (forward + reverse complement)
# For seeds with many Ns, this could be slow since it processes the entire FASTQ for each k-mer pattern sequentially


usage() {
    cat << EOF
Usage: $0 -k <kmers_file> -f <fastq_file> [-o <output_file>] [-h]

Count k-mers with N wildcards in FASTQ files

Arguments:
    -k  K-mers file (one k-mer per line)
    -f  FASTQ file (can be gzipped)
    -o  Output file (default: stdout)
    -h  Show this help message

Example:
    $0 -k kmers.txt -f reads.fastq
    $0 -k kmers.txt -f reads.fastq.gz -o counts.txt
EOF
    exit 1
}

# Parse arguments
while getopts "k:f:o:h" opt; do
    case $opt in
        k) KMERS_FILE="$OPTARG";;
        f) FASTQ_FILE="$OPTARG";;
        o) OUTPUT_FILE="$OPTARG";;
        h) usage;;
        *) usage;;
    esac
done

# Check required arguments
if [ -z "$KMERS_FILE" ] || [ -z "$FASTQ_FILE" ]; then
    echo "Error: Both k-mer file and FASTQ file are required"
    usage
fi

# Check if files exist
if [ ! -f "$KMERS_FILE" ]; then
    echo "Error: K-mer file '$KMERS_FILE' not found"
    exit 1
fi

if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file '$FASTQ_FILE' not found"
    exit 1
fi

# Function to read FASTQ (handles both plain and gzipped)
read_fastq() {
    if [[ "$1" == *.gz ]]; then
        zcat "$1"
    else
        cat "$1"
    fi
}

# Process k-mers
process_kmers() {
    #echo "# K-mer counts for: $FASTQ_FILE"
    #echo "# K-mer"$'\t'"Count"
    
    while read kmer; do
        # Skip empty lines and comments
        [ -z "$kmer" ] && continue
        [[ "$kmer" =~ ^# ]] && continue
        
        # Wildcard handling: 
        # Converts N characters to regex pattern [ATCG], so a k-mer like "ATNCG" becomes "AT[ATCG]CG" for pattern matching.
        # Convert N to [ATCG] for regex
        pattern=$(echo $kmer | sed 's/N/[ATCG]/g')
        
        # Count occurrences
        # - Reads the FASTQ file (every 4th line starting at line 2 contains sequences), Multi-pass: Reads entire FASTQ for EACH k-mer. 
        # O(n × m) complexity where n = number of k-mers, m = FASTQ size. If you have 100 k-mers, it reads the FASTQ file 100 times.
        # - Uses awk to match the pattern against sequence lines
        # - Counts total matches for each k-mer
        count=$(read_fastq "$FASTQ_FILE" | \
                awk -v p="$pattern" 'NR%4==2 && $0 ~ p {c++} END {print c+0}')
        
        echo -e "$kmer\t$count"
    done < "$KMERS_FILE"
}

# Output results
if [ -n "$OUTPUT_FILE" ]; then
    process_kmers > "$OUTPUT_FILE"
    echo "Results saved to: $OUTPUT_FILE"
else
    process_kmers
fi