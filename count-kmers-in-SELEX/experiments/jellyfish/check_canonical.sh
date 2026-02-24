#!/bin/bash

COUNTS_FILE="$1"

# Check if reverse complements are present
check_canonical() {
    # Sample some k-mers and check for their reverse complements
    head -100 "$COUNTS_FILE" | while read kmer count; do
        rc=$(echo "$kmer" | tr 'ATCG' 'TAGC' | rev)
        if [ "$kmer" != "$rc" ]; then
            # Check if RC exists in file
            if grep -q "^$rc " "$COUNTS_FILE"; then
                echo "Found both $kmer and $rc - NOT canonical"
                return 1
            fi
        fi
    done
    echo "Appears to be canonical counts"
    return 0
}

check_canonical