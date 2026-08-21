#!/bin/bash

#matches jellyfish's canonical behavior where each k-mer and its reverse complement are considered equivalent, but palindromes are only counted once.

#./count_kmers_canonical.sh -k $kmers -f $input$sample_base -o $out_path$motif"_signal.txt" 

usage() {
    cat << EOF
Usage: $0 -k <kmers_file> -f <fastq_file> [-o <output_file>] [-h]
Count canonical k-mers with N wildcards in FASTQ files
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
        k) KMERS_FILE="$OPTARG";; #KMERS_FILE=kmers
        f) FASTQ_FILE="$OPTARG";; #FASTQ_FILE=input$sample_base
        o) OUTPUT_FILE="$OPTARG";; #OUTPUT_FILE=test.txt
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

# Function to compute reverse complement
# Reverses the string and complements bases
# N remains N (not included in tr translation)
reverse_complement() {
    echo "$1" | rev | tr 'ATCGatcg' 'TAGCtagc'
}

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
    while read kmer; do
        # Skip empty lines and comments
        [ -z "$kmer" ] && continue
        [[ "$kmer" =~ ^# ]] && continue
        
        # Convert to uppercase
        kmer=$(echo "$kmer" | tr '[:lower:]' '[:upper:]')
        
        # Get reverse complement
        rc_kmer=$(reverse_complement "$kmer")
        
        # Palindrome detection
        # Checks if k-mer equals its reverse complement
        # Example: "ATNTAT" → reverse complement is "ATNTAT" (palindromic)
        # Check if palindromic (self-reverse-complementary)
        if [ "$kmer" = "$rc_kmer" ]; then
            # Palindromic: count only the original
            # Palindromic: Counts only the original pattern
            pattern=$(echo "$kmer" | sed 's/N/[ATCG]/g')
            count=$(read_fastq "$FASTQ_FILE" | \
                    awk -v p="$pattern" 'NR%4==2 && $0 ~ p {c++} END {print c+0}')
        else
            # Non-palindromic: count both forward and reverse complement
            # Non-palindromic: Counts both forward and reverse complement in single AWK pass
            pattern1=$(echo "$kmer" | sed 's/N/[ATCG]/g')
            pattern2=$(echo "$rc_kmer" | sed 's/N/[ATCG]/g')
            
            # Count both patterns in one pass through AWK
            # Checks both patterns in one sequence scan
            # Avoids reading the FASTQ twice per k-mer
            count=$(read_fastq "$FASTQ_FILE" | \
                    awk -v p1="$pattern1" -v p2="$pattern2" '
                    NR%4==2 {
                        if ($0 ~ p1) c1++
                        if ($0 ~ p2) c2++
                    } 
                    END {print c1+c2+0}')
        fi
        
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