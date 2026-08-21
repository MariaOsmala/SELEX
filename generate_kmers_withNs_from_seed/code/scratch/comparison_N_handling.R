#!/usr/bin/env Rscript

# Comparison script: Old vs New N handling
# Shows the difference between treating N as wildcard vs fixed

cat("=======================================================\n")
cat("Comparison: Old vs New handling of N in k-mer generation\n")
cat("=======================================================\n\n")

# Test seed
seed <- "ATNGC"
cat("Test seed:", seed, "\n")
cat("Length:", nchar(seed), "\n")
cat("N at position: 3\n\n")

# OLD BEHAVIOR: N expands to A, C, G, T
cat("OLD BEHAVIOR (N expands to all 4 bases):\n")
cat("-----------------------------------------\n")

old_realizations <- c(
  "ATAGC",  # N -> A
  "ATCGC",  # N -> C
  "ATGGC",  # N -> G
  "ATTGC"   # N -> T
)

cat("Realizations (4 total):\n")
for (r in old_realizations) {
  cat("  ", r, "\n")
}

cat("\nFor each realization, Hamming-1 neighbors would include:\n")
cat("  - Substitutions at positions 1, 2, 3, 4, 5 (all positions)\n")
cat("  - Total: 4 realizations × 5 positions × 3 alternatives = 60 neighbors\n")
cat("  - But position 3 changes in realizations overlap with Hamming-1\n\n")

# NEW BEHAVIOR: N stays as N
cat("NEW BEHAVIOR (N stays as N):\n")
cat("-----------------------------\n")

new_realizations <- c(
  "ATNGC"  # N stays as N
)

cat("Realizations (1 total):\n")
for (r in new_realizations) {
  cat("  ", r, "\n")
}

cat("\nFor this realization, Hamming-1 neighbors:\n")
cat("  - Substitutions ONLY at positions 1, 2, 4, 5 (not position 3 with N)\n")
cat("  - Position 1 (A): Can change to C, G, T\n")
cat("  - Position 2 (T): Can change to A, C, G\n")
cat("  - Position 4 (G): Can change to A, C, T\n")
cat("  - Position 5 (C): Can change to A, G, T\n")
cat("  - Total: 4 positions × 3 alternatives = 12 neighbors\n\n")

cat("Hamming-1 neighbors (12 total):\n")
hamming1_new <- c(
  "CTNGC", "GTNGC", "TTNGC",  # Position 1 changes
  "AANGC", "ACNGC", "AGNGC",  # Position 2 changes
  "ATNAC", "ATNCC", "ATNTC",  # Position 4 changes
  "ATNGA", "ATNGG", "ATNGT"   # Position 5 changes
)
for (h in hamming1_new) {
  cat("  ", h, "\n")
}

cat("\n=======================================================\n")
cat("More complex example with multiple Ns\n")
cat("=======================================================\n\n")

seed2 <- "NNCCGGNNNNNNCCGGNN"
cat("Test seed:", seed2, "\n")
cat("Length:", nchar(seed2), "\n")
cat("N positions: 1, 2, 7-12, 17, 18 (10 Ns total)\n")
cat("Non-N positions: 3-6, 13-16 (8 positions)\n\n")

cat("OLD BEHAVIOR:\n")
cat("  Realizations: 4^10 = 1,048,576 (each N expands to 4 bases)\n")
cat("  Hamming-1: Each realization can mutate at all 18 positions\n\n")

cat("NEW BEHAVIOR:\n")
cat("  Realizations: 1 (all Ns stay as N)\n")
cat("  Single realization: NNCCGGNNNNNNCCGGNN\n")
cat("  Hamming-1: Can only mutate at 8 non-N positions (3-6, 13-16)\n")
cat("  Total Hamming-1: 8 positions × 3 alternatives = 24 neighbors\n\n")

cat("Example Hamming-1 neighbors (first 10):\n")
example_h1 <- c(
  "NNACGGNNNNNNCCGGNN",  # Position 3: C->A
  "NNGCGGNNNNNNCCGGNN",  # Position 3: C->G
  "NNTCGGNNNNNNCCGGNN",  # Position 3: C->T
  "NNCCAGNNNNNNCCGGNN",  # Position 4: G->A
  "NNCCCGNNNNNNCCGGNN",  # Position 4: G->C
  "NNCCGGNNNNNNACGGNN",  # Position 13: C->A
  "NNCCGGNNNNNNCCAGNN",  # Position 14: G->A
  "NNCCGGNNNNNNCCCGNN",  # Position 14: G->C
  "NNCCGGNNNNNNCCGANN",  # Position 15: G->A
  "NNCCGGNNNNNNCCGCNN"   # Position 15: G->C
)
for (i in 1:10) {
  cat("  ", example_h1[i], "\n")
}

cat("\nHamming-2 neighbors:\n")
cat("  Combinations of 2 positions from the 8 mutable positions\n")
cat("  Total combinations: choose(8,2) = 28 pairs\n")
cat("  Each pair has 3×3 = 9 possible double mutations\n")
cat("  Total Hamming-2: 28 × 9 = 252 neighbors\n\n")

cat("=======================================================\n")
cat("KEY DIFFERENCES SUMMARY\n")
cat("=======================================================\n\n")

cat("1. Number of realizations:\n")
cat("   OLD: Exponential in number of Ns (4^n where n = number of Ns)\n")
cat("   NEW: Always 1 (Ns stay fixed)\n\n")

cat("2. Positions available for mutation:\n")
cat("   OLD: All positions can mutate\n")
cat("   NEW: Only non-N positions can mutate\n\n")

cat("3. Computational efficiency:\n")
cat("   OLD: Can be very slow with many Ns (millions of realizations)\n")
cat("   NEW: Fast, as only one realization to process\n\n")

cat("4. Biological relevance:\n")
cat("   NEW: More appropriate when N represents unknown/any base\n")
cat("        that should remain flexible in the final sequences\n")
