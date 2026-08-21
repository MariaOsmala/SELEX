IUPAC_MAP <- Biostrings::IUPAC_CODE_MAP
BASES <- c("A","C","G","T")

.check_iupac <- function(seq) {
  ch <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  bad <- setdiff(unique(ch), names(IUPAC_MAP))
  if (length(bad)) stop("Found non-IUPAC characters: ", paste(bad, collapse=" "))
  ch
}
opts_per_pos <- function(chars) lapply(chars, function(ch) strsplit(IUPAC_MAP[[ch]], "", fixed = TRUE)[[1]])

rev_comp <- function(seq) { paste0(rev(strsplit(chartr("ACGT","TGCA", seq), "", fixed = TRUE)[[1]]), collapse = "") }
canonical_form <- function(seq) { rc <- rev_comp(seq); if (seq <= rc) seq else rc }

expand_iupac_stream <- function(iupac_seq, writer) {
  chars <- .check_iupac(iupac_seq); opts <- opts_per_pos(chars); k <- length(opts); buf <- character(k)
  rec <- function(pos) {
    if (pos > k) { writer(paste0(buf, collapse = "")); return(invisible()) }
    for (b in opts[[pos]]) { buf[pos] <- b; rec(pos + 1L) }
  }
  rec(1L); invisible()
}

emit_hamming1_neighbors <- function(real_seq, allowed_idx, writer) {
  s <- strsplit(real_seq, "", fixed = TRUE)[[1]]
  L <- length(s)
  idx <- as.integer(allowed_idx[!is.na(allowed_idx) & allowed_idx >= 1L & allowed_idx <= L])
  if (!length(idx)) return(0L)
  
  wrote <- 0L
  for (i in idx) {
    old <- s[i]
    if (is.na(old) || !nzchar(old)) next
    for (b in BASES) {
      # NA-proof compare
      if (isTRUE(b != old)) {
        s[i] <- b
        writer(paste0(s, collapse = ""))
        wrote <- wrote + 1L
      }
    }
    s[i] <- old
  }
  wrote
}



process_seed <- function(seed,
                         out_realizations = NULL,
                         out_hamming1     = NULL,
                         compress = TRUE,
                         canonicalize = TRUE,
                         distinct_realizations = FALSE,
                         distinct_hamming1     = FALSE,
                         log_first_h1 = 3) {
  
  seedU <- toupper(seed)
  chars <- .check_iupac(seedU)
  k <- length(chars)
  
  # positions where the ORIGINAL seed is NOT 'N'
  allowed_idx <- which(chars != "N")
  allowed_idx <- as.integer(allowed_idx[!is.na(allowed_idx) & allowed_idx >= 1L & allowed_idx <= k])
  
  message("Mutation mask (non-N positions): ", if (length(allowed_idx)) paste(allowed_idx, collapse = ",") else "<none>")
  
  # fail fast if there are no mutable positions — prevents H1=0 surprises
  if (!length(allowed_idx)) {
    message("No non-N positions in the seed -> no Hamming-1 neighbors by design.")
  }
  
  conR <- conH <- NULL
  on.exit({ if (!is.null(conR)) close(conR); if (!is.null(conH)) close(conH) }, add = TRUE)
  
  if (!is.null(out_realizations)) {
    conR <- if (compress) gzfile(paste0(out_realizations, ".gz"), "wt") else file(out_realizations, "wt")
  }
  if (!is.null(out_hamming1)) {
    conH <- if (compress) gzfile(paste0(out_hamming1, ".gz"), "wt") else file(out_hamming1, "wt")
  }
  
  seenR <- if (distinct_realizations) new.env(hash = TRUE, parent = emptyenv(), size = 1e6) else NULL
  seenH <- if (distinct_hamming1)     new.env(hash = TRUE, parent = emptyenv(), size = 1e6) else NULL
  
  writeR <- if (is.null(conR)) function(x){} else function(x) writeLines(x, con = conR, useBytes = TRUE)
  writeH <- if (is.null(conH)) function(x){} else function(x) writeLines(x, con = conH, useBytes = TRUE)
  
  .norm_key <- function(x) { k <- as.character(x)[1L]; if (!is.character(k) || is.na(k) || !nzchar(k)) NA_character_ else k }
  
  real_count <- 0L
  h1_attempt <- 0L  # before canonical/dedup
  h1_written <- 0L  # after canonical/dedup
  
  canon_writer_R <- function(seq) {
    s <- if (canonicalize) canonical_form(seq) else seq
    if (!is.null(seenR)) { key <- .norm_key(s); if (!is.na(key)) { if (!is.null(seenR[[key]])) return(invisible(NULL)); seenR[[key]] <- TRUE } }
    writeR(s); real_count <<- real_count + 1L
  }
  
  
  # --- Replace your canon_writer_H with this (no local(), clear counters) ---
  # Put this definition INSIDE process_seed() after you declare
  #   h1_attempt <- 0L; h1_written <- 0L
  canon_writer_H <- {
    first_dump <- 0L
    function(seq) {
      h1_attempt <<- h1_attempt + 1L
      s <- if (canonicalize) canonical_form(seq) else seq
      if (!is.null(seenH)) {
        key <- .norm_key(s)
        if (!is.na(key) && !is.null(seenH[[key]])) return(invisible(NULL))
        if (!is.na(key)) seenH[[key]] <- TRUE
      }
      writeH(s)
      h1_written <<- h1_written + 1L
      if (log_first_h1 > 0 && first_dump < log_first_h1) {
        message("H1 example: ", s)
        first_dump <<- first_dump + 1L
      }
      invisible(NULL)
    }
  }
  
  
  # ALWAYS call the H1 emitter; writer itself is a no-op if no H1 file is open
  expand_iupac_stream(seedU, function(real) {
    canon_writer_R(real)
    emit_hamming1_neighbors(real, allowed_idx, canon_writer_H)
  })
  
  message(sprintf("Done: realizations=%s | hamming1_attempted=%s | hamming1_written=%s",
                  format(real_count, big.mark=","), format(h1_attempt, big.mark=","), format(h1_written, big.mark=",")))
  invisible(list(realizations = real_count, h1_attempted = h1_attempt, h1_written = h1_written,
                 mutable_positions = allowed_idx))
}
    

# ---------- Example ----------
library("dplyr")
library("readr")
metadata=read_tsv("/projappl/project_2013895/SELEX/generate_kmers_from_seed/code/metadata_failed_and_successfull_seeds_degeneracy.tsv")


motif=metadata %>% filter(kmer_generation_successfull==TRUE) %>% filter(length==(min(length)+1)) %>% pull(ID)
motif=motif[1]
s=metadata %>% filter(ID==motif) %>% pull(seed)

print(seq_stats(s))
# $ k                : int 7
# $ mutable_positions: int [1:5] 2 3 4 5 6
# $ n_realizations   : num 16 4^2
# $ hamming1_per_real: int 15 3*5
# $ total_hamming1   : num 240 =16*15

#s=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(seed)
#motif=metadata %>% filter(kmer_generation_successfull==FALSE) %>% filter(length==min(length)) %>% pull(ID)
#print(seq_stats(s))
 
 #List of 5
 #$ k                : int 18
 #$ mutable_positions: int [1:8] 3 4 5 6 13 14 15 16
 #$ n_realizations   : num 1048576
 #$ hamming1_per_real: int 24
 #$ total_hamming1   : num 25165824
 
 
process_seed(
  seed=s,
  out_realizations = paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif,"_realizations.txt"),
  out_hamming1     = paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif,"_h1.txt"),
  compress = TRUE,
  canonicalize = TRUE,
  distinct_realizations = TRUE,  # set TRUE to dedup canonicals in-memory
  distinct_hamming1     = TRUE
)


process_seed("AGATAANN",
             out_realizations = paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif,"_realizations.txt"),
             out_hamming1     = paste0("/scratch/project_2013895/SELEX/kmers_from_long_degenerate_seeds/", motif,"_h1.txt"),
             compress = FALSE,
             canonicalize = TRUE,
             distinct_realizations = FALSE,
             distinct_hamming1 = FALSE)

# seeds <- c("NNCCGGNNNNNNCCGGNN", "NRYGKMBDHVN")
# process_seeds(seeds, output_dir = "out/", prefix = "kmers",
#               compress = TRUE, canonicalize = TRUE)
