# Function to align one sequence against another (considering reverse complement)
align_with_revcomp <- function(seq1, seq2, substitutionMatrix = NULL) {
  if (is.null(substitutionMatrix)) {
    substitutionMatrix <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  }
  
  # Forward alignment
  aln_forward <- pwalign::pairwiseAlignment(seq1, seq2, substitutionMatrix = substitutionMatrix,
                                            type="local", #Smith-Waterman
                                            gapOpening = -2, gapExtension = -1)
  
  # Reverse complement alignment
  aln_rev <- pwalign::pairwiseAlignment(seq1, reverseComplement(seq2), substitutionMatrix = substitutionMatrix,
                                        type="local", #Smith-Waterman
                                        gapOpening = -2, gapExtension = -1)
  
  # Return whichever has higher score
  if (score(aln_rev) > score(aln_forward)) {
    return(list(best = aln_rev, forward = aln_forward, reverse = aln_rev))
  } else {
    return(list(best = aln_forward, forward = aln_forward, reverse = aln_rev))
  }
}
