# BiocManager::install(c("Biostrings","ShortRead","BiocParallel"))
library(Biostrings); library(ShortRead); library(BiocParallel)

iupac_count <- function(patterns, fastqs, merge_strands=FALSE, chunk=1e6, BPPARAM=MulticoreParam(8)) {
  pats <- DNAStringSet(patterns); rc <- reverseComplement(pats)
  is_pal <- as.character(pats) == as.character(rc)
  totals <- integer(length(pats)); names(totals) <- patterns
  
  for (fq in fastqs) {
    #fq=fastqs[1]
    st <- FastqStreamer(fq, n=chunk); on.exit(close(st), add=TRUE)
    repeat {
      fqch <- yield(st); if (length(fqch)==0) break
      reads <- sread(fqch)
      # parallel over patterns
      counts <- bplapply(seq_along(pats), function(i){
        c1 <- sum(vcountPattern(pats[i], reads, fixed=FALSE, with.indels=FALSE))
        if (!merge_strands) return(c1)
        if (is_pal[i]) return(c1)
        c2 <- sum(vcountPattern(rc[i], reads, fixed=FALSE, with.indels=FALSE))
        c1 + c2
      }, BPPARAM=BPPARAM)
      totals <- totals + unlist(counts)
    }
  }
  data.frame(pattern=names(totals), count=totals)
}

patterns=read.table("/scratch/project_2013895/SELEX/combined_from_fixed_N_seeds/Meis2_HT-SELEX_TCAAAA20NTA_P_NTGACAN_1_3.txt")
patterns=XString(patterns$V1)
fastqs="/scratch/project_2013895/SELEX/data_exclude_reads_with_N/Meis2_TCAAAA20NTA_P_3.fastq.gz"
