## =======================
## Utilities: IUPAC & Hamming neighbors
## =======================

library("tidyverse")
library("dplyr")

.libPaths("/projappl/project_2013895/softwares/RPackages")
#install.packages("RcppAlgos")
library("RcppAlgos")

# install.packages("digest")  # if needed
#library(digest)

## =======================
## Parsing your “ID ↔ seed” list
## =======================

# Option A: if the content is already in a character vector `raw_lines`
# raw_lines <- readLines("your_file.txt")

# Using the pasted content directly: replace the next line with readLines(...) in practice

args <- commandArgs(trailingOnly = TRUE)
arrays <- as.numeric(args[1]) # 0-363


start_ind=arrays*10+1 #100
end_ind=(arrays+1)*10 #100


#3636 motifs
if(end_ind>3636){ # 
  end_ind=3636
}


metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

metadata=metadata %>% filter(!(experiment %in% "Methyl-HT-SELEX")) #3636



seeds <- toupper(metadata$seed)
names(seeds)=metadata$ID





.bitset_create <- function(nbits) {
  nbytes <- as.integer(ceiling(nbits / 8))
  structure(raw(nbytes), class = "bitset", nbits = as.double(nbits))
}
.bitset_set <- function(bs, i) {
  byte <- (i %/% 8L) + 1L
  off  <- i %% 8L
  bs[byte] <- as.raw(bitwOr(as.integer(bs[byte]), bitwShiftL(1L, off)))
  bs
}
.bitset_test <- function(bs, i) {
  byte <- (i %/% 8L) + 1L
  off  <- i %% 8L
  bitwAnd(as.integer(bs[byte]), bitwShiftL(1L, off)) != 0L
}

# ---- Bloom filter ----
bloom_create <- function(n_expected, fp = 1e-3) {
  m_bits <- ceiling(-(n_expected * log(fp)) / (log(2)^2))
  k_hash <- max(1L, round((m_bits / n_expected) * log(2)))
  list(
    m = as.integer(m_bits),
    k = as.integer(k_hash),
    bits = .bitset_create(m_bits)
  )
}

# Compute (digest(key) mod m) robustly from raw bytes (no overflow)
.digest_mod <- function(key, salt, m) {
  # xxhash32 as raw (4 bytes)
  r <- digest::digest(paste0(salt, key), algo = "xxhash32", serialize = FALSE, raw = TRUE)
  # Horner: ((((b0*256 + b1)*256 + b2)*256 + b3) %% m)
  acc <- 0L
  for (b in as.integer(r)) {
    acc <- (acc * 256L + b) %% m
  }
  acc
}

# Double hashing: indices = (h1 + j*h2) %% m
.bloom_indices <- function(key, m, k) {
  h1 <- .digest_mod(key, "1", m)
  h2 <- .digest_mod(key, "2", m)
  if (h2 == 0L) h2 <- 1L
  # build k indices
  (h1 + (0:(k-1L)) * h2) %% m
}

bloom_maybe_add <- function(bf, key) {
  idxs <- .bloom_indices(key, bf$m, bf$k)
  # all indices are integers in [0, m); no NA
  all_set <- TRUE
  for (idx in idxs) if (!.bitset_test(bf$bits, idx)) { all_set <- FALSE; break }
  if (all_set) return(FALSE)   # probably seen
  for (idx in idxs) bf$bits <<- .bitset_set(bf$bits, idx)
  TRUE                          # probably new
}


.get_iupac_map <- function() {
  # Prefer Biostrings’ authoritative map if available
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    m <- Biostrings::IUPAC_CODE_MAP
    return(lapply(m, function(s) strsplit(s, "", fixed = TRUE)[[1]]))
  }
  # Fallback (DNA only)
  fallback <- c(
    A="A", C="C", G="G", T="T",
    R="AG", Y="CT", S="GC", W="AT", K="GT", M="AC",
    B="CGT", D="AGT", H="ACT", V="ACG",
    N="ACGT"
  )
  lapply(fallback, function(s) strsplit(s, "", fixed = TRUE)[[1]])
}


choices_from_seed <- function(seed, alphabet = c("A","C","G","T")) {
  mp <- .get_iupac_map()
  chars <- strsplit(toupper(seed), "", fixed = TRUE)[[1]]
  lapply(chars, function(ch) {
    if (!ch %in% names(mp)) stop(sprintf("Unknown IUPAC symbol: %s", ch))
    x <- intersect(mp[[ch]], alphabet)
    if (!length(x)) stop(sprintf("Symbol %s incompatible with alphabet", ch))
    x
  })
}


# Calls `handler(batch)` repeatedly with batches of concrete expansions.
expand_iupac_stream <- function(seed, alphabet = c("A","C","G","T"),
                                handler, batch_size = 100000L) {
  choices <- choices_from_seed(seed, alphabet)
  k <- length(choices)
  lens <- vapply(choices, length, integer(1))
  # Early exit if any position had 0 options
  if (any(lens == 0L)) return(invisible(NULL))
  
  idx <- rep(1L, k)                 # 1-based indices into each choices[[i]]
  buf <- character(0L)
  letters_here <- character(k)
  
  repeat {
    # Build one string from the current index tuple
    for (i in seq_len(k)) letters_here[i] <- choices[[i]][idx[i]]
    buf <- c(buf, paste0(letters_here, collapse = ""))
    
    # Flush a full batch
    if (length(buf) >= batch_size) { handler(buf); buf <- character(0L) }
    
    # Odometer increment (rightmost position carries left)
    i <- k
    while (i >= 1L && idx[i] == lens[i]) { idx[i] <- 1L; i <- i - 1L }
    if (i < 1L) break  # done
    idx[i] <- idx[i] + 1L
  }
  
  if (length(buf)) handler(buf)      # flush tail
  invisible(NULL)
}




expand_iupac <- function(seed, alphabet = c("A","C","G","T")) {
  seed <- toupper(seed)
  mp <- .get_iupac_map()
  chars <- strsplit(seed, "", fixed = TRUE)[[1]]
  choices <- lapply(chars, function(ch) {
    if (!ch %in% names(mp)) stop(sprintf("Unknown IUPAC symbol: %s", ch))
    x <- intersect(mp[[ch]], alphabet)
    if (!length(x)) stop(sprintf("Symbol %s incompatible with alphabet", ch))
    x
  })
  # Fast Cartesian product; switch to RcppAlgos if available for large expansions
  if (requireNamespace("RcppAlgos", quietly = TRUE)) {
    grid <- do.call(RcppAlgos::expandGrid, choices)  # matrix (rows = products)
    apply(grid, 1, paste0, collapse = "")
  } else {
    res <- ""  # iterative outer-product in base R
    for (opts in choices) res <- as.vector(outer(res, opts, paste0))
    res
  }
}


# Generic product streamer over a list of per-position choices
product_stream <- function(choices, handler, batch_size = 100000L) {
  k <- length(choices)
  lens <- vapply(choices, length, integer(1))
  if (any(lens == 0L)) return(invisible(NULL))
  idx <- rep(1L, k)
  letters_here <- character(k)
  buf <- character(0L)
  repeat {
    for (i in seq_len(k)) letters_here[i] <- choices[[i]][idx[i]]
    buf <- c(buf, paste0(letters_here, collapse = ""))
    if (length(buf) >= batch_size) { handler(buf); buf <- character(0L) }
    i <- k
    while (i >= 1L && idx[i] == lens[i]) { idx[i] <- 1L; i <- i - 1L }
    if (i < 1L) break
    idx[i] <- idx[i] + 1L
  }
  if (length(buf)) handler(buf)
  invisible(NULL)
}



#stream the seed’s concrete realizations (d=0)
stream_deg_expansion <- function(seed, alphabet = c("A","C","G","T"),
                                 handler, batch_size = 100000L) {
  A <- choices_from_seed(seed, alphabet)
  product_stream(A, handler, batch_size)
}

stream_deg_hamming1_from_choices <- function(A, handler, alphabet = c("A","C","G","T"),
                                             batch_size = 100000L) {
  k <- length(A)
  for (i in seq_len(k)) {
    Ci <- setdiff(alphabet, A[[i]]); if (!length(Ci)) next
    B <- A; B[[i]] <- Ci                          # outside at i, inside elsewhere
    product_stream(B, handler, batch_size)
  }
}

stream_deg_hamming2_from_choices <- function(A, handler, alphabet = c("A","C","G","T"),
                                             batch_size = 100000L) {
  k <- length(A); if (k < 2) return(invisible(NULL))
  pr <- utils::combn(k, 2L)
  for (c in seq_len(ncol(pr))) {
    i <- pr[1,c]; j <- pr[2,c]
    Ci <- setdiff(alphabet, A[[i]]); if (!length(Ci)) next
    Cj <- setdiff(alphabet, A[[j]]); if (!length(Cj)) next
    B <- A; B[[i]] <- Ci; B[[j]] <- Cj            # outside at i & j, inside elsewhere
    product_stream(B, handler, batch_size)
  }
}

neighbors_hamming1 <- function(kmer, alphabet = c("A","C","G","T")) {
  s <- strsplit(kmer, "", fixed = TRUE)[[1]]
  k <- length(s)
  out <- vector("list", k)
  for (i in seq_len(k)) {
    repl <- setdiff(alphabet, s[i])
    if (length(repl)) {
      prefix <- if (i > 1) paste(s[1:(i-1)], collapse = "") else ""
      suffix <- if (i < k) paste(s[(i+1):k], collapse = "") else ""
      out[[i]] <- paste0(prefix, repl, suffix)
    } else out[[i]] <- character()
  }
  unique(unlist(out, use.names = FALSE))
}

neighbors_hamming2 <- function(kmer, alphabet = c("A","C","G","T")) {
  s <- strsplit(kmer, "", fixed = TRUE)[[1]]
  k <- length(s)
  if (k < 2) return(character())
  pairs <- utils::combn(k, 2)
  out <- vector("list", ncol(pairs))
  for (j in seq_len(ncol(pairs))) {
    i <- pairs[1, j]; h <- pairs[2, j]
    prefix    <- if (i > 1) paste(s[1:(i-1)], collapse = "") else ""
    mid_block <- if (h > i+1) paste(s[(i+1):(h-1)], collapse = "") else ""
    suffix    <- if (h < k) paste(s[(h+1):k], collapse = "") else ""
    repl_i <- setdiff(alphabet, s[i])
    repl_h <- setdiff(alphabet, s[h])
    if (length(repl_i) && length(repl_h)) {
      out[[j]] <- as.vector(outer(repl_i, repl_h, \(a,b) paste0(prefix, a, mid_block, b, suffix)))
    } else out[[j]] <- character()
  }
  unique(unlist(out, use.names = FALSE))
}

kmers_from_seed_with_neighbors <- function(seed,
                                           alphabet = c("A","C","G","T"),
                                           include_seed_expansion = TRUE,
                                           include_d1 = TRUE,
                                           include_d2 = TRUE,
                                           dedup = TRUE) {
  expanded <- expand_iupac(seed, alphabet)
  res <- character()
  if (include_seed_expansion) res <- c(res, expanded)
  if (include_d1) for (s in expanded) res <- c(res, neighbors_hamming1(s, alphabet))
  if (include_d2) for (s in expanded) res <- c(res, neighbors_hamming2(s, alphabet))
  if (dedup) res <- unique(res)
  res
}

kmers_from_seed_with_neighbors_canonical <- function(seed,
                                                     alphabet = c("A","C","G","T"),
                                                     include_seed_expansion = TRUE,
                                                     include_d1 = TRUE,
                                                     include_d2 = TRUE, dedup=FALSE) {
  # generate (as before)
  km <- kmers_from_seed_with_neighbors(seed, alphabet,
                                       include_seed_expansion, include_d1, include_d2,
                                       dedup = FALSE)  # defer dedup
  unique(canonicalize(km))
}



## Reverse complement (vectorized over character vector of A/C/G/T strings)
rc <- function(x) {
  # complement using chartr, then reverse each string
  comp <- chartr("ACGT", "TGCA", x)
  vapply(strsplit(comp, "", fixed = TRUE), function(cs) paste0(rev(cs), collapse = ""), "")
}


## Canonicalize: lexicographic min of s and rc(s)
canonicalize <- function(v) {
  rcv <- rc(v)
  # pmin is vectorized lexicographic min for same-length strings
  pmin(v, rcv)
}

## Write a character vector to a connection in chunks
.write_lines_chunked <- function(v, con, chunk = 500000L) {
  n <- length(v); if (n == 0L) return(invisible())
  i <- 1L
  while (i <= n) {
    j <- min(i + chunk - 1L, n)
    writeLines(v[i:j], con, sep = "\n", useBytes = TRUE)
    i <- j + 1L
  }
}




## =======================
## Run generation
## =======================

# Choose output mode:
#per_seed_results <- TRUE      # keep neighbors per seed separately
#write_per_seed   <- TRUE     # set TRUE to stream each seed's results to disk (fast & RAM-safe)

# Directory for streaming output (only used if write_per_seed = TRUE)
results_path <- "/scratch/project_2013895/SELEX/seed_kmers_out"

out_dir <- Sys.getenv("LOCAL_SCRATCH")

#if (write_per_seed && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alphabet <- c("A","C","G","T")

# Function to generate (streaming optional)
generate_for_one_seed <- function(seed_str, seed_id) {
  if (write_per_seed) {
    fn <- file.path(out_dir, paste0(seed_id, "_kmers_d0_d1_d2.txt.gz"))
    con <- gzfile(fn, open = "w")
    on.exit(close(con), add = TRUE)
    # Stream: write expansion, then d1, then d2 (dedup at the end if you like via sort|uniq externally)
    expanded <- expand_iupac(seed_str, alphabet)
    writeLines(expanded, con)
    for (s in expanded) writeLines(neighbors_hamming1(s, alphabet), con)
    for (s in expanded) writeLines(neighbors_hamming2(s, alphabet), con)
    return(invisible(fn))
  } else {
    km <- kmers_from_seed_with_neighbors_canonical(seed_str, alphabet,
                                         include_seed_expansion = TRUE,
                                         include_d1 = TRUE,
                                         include_d2 = TRUE,
                                         dedup = TRUE)
    km
  }
}

make_emitter <- function(con) {
  if (dedup_canonical_in_memory) {
    seen <- new.env(parent = emptyenv())
    function(x) {
      if (!length(x)) return()
      x <- canonicalize(x)
      keep <- !vapply(x, exists, logical(1), envir = seen, USE.NAMES = FALSE)
      if (any(keep)) {
        invisible(lapply(x[keep], function(k) assign(k, TRUE, envir = seen)))
        .write_lines_chunked(x[keep], con)
      }
    }
  } else {
    function(x) { if (length(x)) .write_lines_chunked(canonicalize(x), con) }
  }
}

write_deg_hamming012 <- function(seed, seed_id,
                                 out_dir,
                                 alphabet = c("A","C","G","T"),
                                 batch_size = 100000L) {
  fn <- file.path(out_dir, paste0(seed_id, "_deg_hamming012.canonical.txt.gz"))
  con <- gzfile(fn, "w"); on.exit(close(con), add = TRUE)
  emit <- make_emitter(con)
  
  # d = 0 (optional)
  stream_deg_expansion(seed, alphabet, emit, batch_size)
  # d = 1
  stream_deg_hamming1(seed, alphabet, emit, batch_size)
  # d = 2
  stream_deg_hamming2(seed, alphabet, emit, batch_size)
  
  invisible(fn)
}

degenerate_counts <- function(seed, sigma = 4L, alphabet = c("A","C","G","T")) {
  A <- choices_from_seed(seed, alphabet)
  m <- vapply(A, length, integer(1))
  E  <- prod(as.double(m))
  N1 <- E * sum((sigma - m) / m)
  if (length(m) >= 2) {
    idx <- utils::combn(seq_along(m), 2L)
    N2 <- sum(E * ((sigma - m[idx[1,]]) * (sigma - m[idx[2,]]) /
                     (m[idx[1,]] * m[idx[2,]])))
  } else N2 <- 0
  list(E = E, N1 = N1, N2 = N2)
}

# toggles
canonicalize_on_write <- TRUE      # << ensure RC-collapsed output
dedup_in_memory       <- FALSE     # TRUE = keep a hash set per seed (uses RAM)
# FALSE = allow dups in file; later do: zcat | sort -u | bgzip
chunk_size            <- 200000L   # adjust for your RAM/throughput

generate_for_one_seed_better <- function(seed_str, seed_id) {
  #seed_str=seeds[4]
  #seed_id=names(seeds)[4]
  fn <- file.path(out_dir, paste0(seed_id, "_kmers_d0_d1_d2.txt.gz"))
  con <- gzfile(fn, open = "w")
  on.exit(close(con), add = TRUE)
  
  # Optional per-seed hash set for de-dup (canonical keys)
  if (dedup_in_memory) seen <- new.env(parent = emptyenv())
  
  emit <- function(x) {
    if (!length(x)) return()
    if (canonicalize_on_write) x <- canonicalize(x)
    if (dedup_in_memory) {
      keep <- !vapply(x, exists, logical(1), envir = seen, USE.NAMES = FALSE)
      if (any(keep)) {
        # mark as seen
        invisible(lapply(x[keep], function(k) assign(k, TRUE, envir = seen)))
        .write_lines_chunked(x[keep], con, chunk = chunk_size)
      }
    } else {
      .write_lines_chunked(x, con, chunk = chunk_size)
    }
  }
  
  # Expand IUPAC seed (vectorized)
  expanded <- expand_iupac(seed_str, alphabet)
  emit(expanded)
  
  # Hamming-1 (stream per expanded k-mer to limit RAM)
  for (s in expanded) emit(neighbors_hamming1(s, alphabet))
  
  # Hamming-2
  for (s in expanded) emit(neighbors_hamming2(s, alphabet))
  system(paste("cp", shQuote(fn), shQuote(results_path)))
  invisible(fn)
}


generate_for_one_seed_better_stream <- function(seed_str, seed_id) {
  #seed_str=seeds[4]
  #seed_id=names(seeds)[4]
  fn <- file.path(out_dir, paste0(seed_id, "_kmers_d0_d1_d2.txt.gz"))
  con <- gzfile(fn, open = "w")
  on.exit(close(con), add = TRUE)
  
  # Optional per-seed hash set for de-dup (canonical keys)
  if (dedup_in_memory) seen <- new.env(parent = emptyenv())
  
  emit <- function(x) {
    if (!length(x)) return()
    if (canonicalize_on_write) x <- canonicalize(x)
    if (dedup_in_memory) {
      keep <- !vapply(x, exists, logical(1), envir = seen, USE.NAMES = FALSE)
      if (any(keep)) {
        # mark as seen
        invisible(lapply(x[keep], function(k) assign(k, TRUE, envir = seen)))
        .write_lines_chunked(x[keep], con, chunk = chunk_size)
      }
    } else {
      .write_lines_chunked(x, con, chunk = chunk_size)
    }
  }
  
  # Stream expansions in batches; for each batch, also emit Hamming-1 & Hamming-2
  expand_iupac_stream(seed_str, alphabet,
                      handler = function(batch) {
                        emit(batch)
                        # Hamming-1
                        for (s in batch) emit(neighbors_hamming1(s, alphabet))
                        # Hamming-2
                        for (s in batch) emit(neighbors_hamming2(s, alphabet))
                      },
                      batch_size = 100000L
  )
  
  system(paste("cp", shQuote(fn), shQuote(results_path)))
  invisible(fn)
}


generate_for_one_seed_direct<- function(seed_str, seed_id,
                                         include_d0 = TRUE, include_d1 = TRUE, include_d2 = TRUE) {
  fn <- file.path(out_dir, paste0(seed_id, "_kmers_d",
                                  paste(which(c(include_d0, include_d1, include_d2)) - 1, collapse = ""),
                                  ".txt.gz"))
  con <- gzfile(fn, open = "w")
  on.exit(close(con), add = TRUE)
  
  # exact per-seed dedup of canonicals (optional)
  if (dedup_in_memory) seen <- new.env(parent = emptyenv())
  
  emit <- function(x) {
    if (!length(x)) return()
    if (canonicalize_on_write) x <- canonicalize(x)
    if (dedup_in_memory) {
      keep <- !vapply(x, exists, logical(1), envir = seen, USE.NAMES = FALSE)
      if (any(keep)) {
        invisible(lapply(x[keep], function(k) assign(k, TRUE, envir = seen)))
        .write_lines_chunked(x[keep], con, chunk = chunk_size)
      }
    } else {
      .write_lines_chunked(x, con, chunk = chunk_size)
    }
  }
  
  A <- choices_from_seed(seed_str, alphabet)
  
  if (include_d0) product_stream(A, emit, batch_size = 100000L)                          # exact d=0
  if (include_d1) stream_deg_hamming1_from_choices(A, emit, alphabet, batch_size = 100000L)  # exact d=1
  if (include_d2) stream_deg_hamming2_from_choices(A, emit, alphabet, batch_size = 100000L)  # exact d=2
  
  # (optional) copy result elsewhere
  if (exists("results_path") && nzchar(results_path)) {
    system(paste("cp", shQuote(fn), shQuote(results_path)))
  }
  invisible(fn)
}



# ---- PARAMETERS you can tweak ----
# alphabet <- c("A","C","G","T")
# bloom_fp          <- 1e-3        # target false positive rate
# bloom_n_expected  <- 5e7         # expected # of canonical lines per seed file (set per your scale)
# chunk_size        <- 200000L     # write chunk size

generate_for_one_seed_bloom <- function(seed_str, seed_id) {
  #i=1
  #seed_str=seeds[i]
  #seed_id=names(seeds)[i]
  
  fn <- file.path(out_dir, paste0(seed_id, "_kmers_d0_d1_d2.canonical.bloom.txt.gz"))
  con <- gzfile(fn, open = "w")
  on.exit(close(con), add = TRUE)
  
  # Bloom filter for this seed
  bf <- bloom_create(n_expected = bloom_n_expected, fp = bloom_fp)
  
  emit <- function(x) {
    if (!length(x)) return()
    x <- canonicalize(x)
    # filter through Bloom: keep only "probably new"
    keep <- vapply(x, function(k) bloom_maybe_add(bf, k), logical(1), USE.NAMES = FALSE)
    if (any(keep)) .write_lines_chunked(x[keep], con, chunk = chunk_size)
  }
  
  # 1) seed expansion (vectorized)
  expanded <- expand_iupac(seed_str, alphabet)
  emit(expanded)
  
  # 2) Hamming-1
  for (s in expanded) emit(neighbors_hamming1(s, alphabet))
  
  # 3) Hamming-2
  for (s in expanded) emit(neighbors_hamming2(s, alphabet))
  
  invisible(fn)
}





seeds=seeds[start_ind:end_ind]

# Run
# if (per_seed_results && !write_per_seed) {
#   results <- lapply(seq_along(seeds), function(i) generate_for_one_seed(seeds[i], names(seeds)[i]))
#   names(results) <- names(seeds)
# } else if (write_per_seed) {
#   invisible(lapply(seq_along(seeds), function(i) generate_for_one_seed(seeds[i], names(seeds)[i])))
#   results <- NULL
# }


dedup_in_memory <- FALSE
invisible(lapply(seq_along(seeds), function(i) generate_for_one_seed_better(seeds[i], names(seeds)[i])))

generate_for_one_seed_better(seeds[1], names(seeds)[1])
generate_for_one_seed_better(seeds[2], names(seeds)[2])
generate_for_one_seed_better(seeds[3], names(seeds)[3])
generate_for_one_seed_better(seeds[4], names(seeds)[4])
generate_for_one_seed_better(seeds[5], names(seeds)[5])
generate_for_one_seed_better(seeds[6], names(seeds)[6])
generate_for_one_seed_better(seeds[7], names(seeds)[7])
generate_for_one_seed_better(seeds[8], names(seeds)[8])
generate_for_one_seed_better(seeds[9], names(seeds)[9])
generate_for_one_seed_better(seeds[10], names(seeds)[10])


#zcat 1234_kmers_d0_d1_d2.txt.gz | sort -u | bgzip > 1234_kmers_canonical_uniq.txt.gz

# sort -u BCL6B_HT-SELEX_TGCGGG20NGA_AC_TGCTTTCTAGGAATTMM_2_4_kmers_d0_d1_d2.txt > BCL6B_HT-SELEX_TGCGGG20NGA_AC_TGCTTTCTAGGAATTMM_2_4_kmers_d0_d1_d2_uniq.txt
# 
# sort -u EGR2_HT-SELEX_TCGGCC20NGA_W_MCGCCCACGCA_2_3_kmers_d0_d1_d2.txt >EGR2_HT-SELEX_TCGGCC20NGA_W_MCGCCCACGCA_2_3_kmers_d0_d1_d2_uniq.txt
# sort -u CTCF_HT-SELEX_TAGCGA20NGCT_AJ_NGCGCCMYCTAGYGGTN_2_4_kmers_d0_d1_d2.txt > CTCF_HT-SELEX_TAGCGA20NGCT_AJ_NGCGCCMYCTAGYGGTN_2_4_kmers_d0_d1_d2_uniq.txt
# sort -u Egr3_HT-SELEX_TCATGC20NGTA_Z_NACGCCCACGCANNN_2_3_kmers_d0_d1_d2.txt > Egr3_HT-SELEX_TCATGC20NGTA_Z_NACGCCCACGCANNN_2_3_kmers_d0_d1_d2_uniq.txt
# sort -u EGR1_HT-SELEX_TACTAT20NATC_AA_NACGCCCACGCANN_2_4_kmers_d0_d1_d2.txt > EGR1_HT-SELEX_TACTAT20NATC_AA_NACGCCCACGCANN_2_4_kmers_d0_d1_d2_uniq.txt
# sort -u EGR3_HT-SELEX_TTAAGA20NGAA_AF_NNACGCCCACGCANN_2_3_kmers_d0_d1_d2.txt > EGR3_HT-SELEX_TTAAGA20NGAA_AF_NNACGCCCACGCANN_2_3_kmers_d0_d1_d2_uniq.txt
# sort -u EGR1_HT-SELEX_TCTCTT20NGA_Y_NMCGCCCMCGCANN_2_2_kmers_d0_d1_d2.txt > EGR1_HT-SELEX_TCTCTT20NGA_Y_NMCGCCCMCGCANN_2_2_kmers_d0_d1_d2_uniq.txt
# sort -u EGR4_HT-SELEX_TAGCGA40NGCT_AI_NNMCGCCCACGCANNN_2_4_kmers_d0_d1_d2.txt > EGR4_HT-SELEX_TAGCGA40NGCT_AI_NNMCGCCCACGCANNN_2_4_kmers_d0_d1_d2_uniq.txt
# sort -u EGR2_HT-SELEX_TCATCA20NCTA_W_NNMCGCCCACGCANN_2_2_kmers_d0_d1_d2.txt > EGR2_HT-SELEX_TCATCA20NCTA_W_NNMCGCCCACGCANN_2_2_kmers_d0_d1_d2_uniq.txt
# sort -u GLI2_HT-SELEX_TAAGTA40NTGA_AI_GACCACMCACNNNG_1_3_kmers_d0_d1_d2.txt > GLI2_HT-SELEX_TAAGTA40NTGA_AI_GACCACMCACNNNG_1_3_kmers_d0_d1_d2_uniq.txt

# 4336 BCL6B_HT-SELEX_TGCGGG20NGA_AC_TGCTTTCTAGGAATTMM_2_4_score_rank.tsv
# 95360 CTCF_HT-SELEX_TAGCGA20NGCT_AJ_NGCGCCMYCTAGYGGTN_2_4_score_rank.tsv
# 33856 EGR1_HT-SELEX_TACTAT20NATC_AA_NACGCCCACGCANN_2_4_score_rank.tsv
# 104704 EGR1_HT-SELEX_TCTCTT20NGA_Y_NMCGCCCMCGCANN_2_2_score_rank.tsv
# 239104 EGR2_HT-SELEX_TCATCA20NCTA_W_NNMCGCCCACGCANN_2_2_score_rank.tsv
# 934 EGR2_HT-SELEX_TCGGCC20NGA_W_MCGCCCACGCA_2_3_score_rank.tsv
# 135424 Egr3_HT-SELEX_TCATGC20NGTA_Z_NACGCCCACGCANNN_2_3_score_rank.tsv
# 135424 EGR3_HT-SELEX_TTAAGA20NGAA_AF_NNACGCCCACGCANN_2_3_score_rank.tsv
# 59776 GLI2_HT-SELEX_TAAGTA40NTGA_AI_GACCACMCACNNNG_1_3_score_rank.tsv


# tmp=sapply(results, length)
# names(tmp)=names(seeds)
# i=1
# result=generate_for_one_seed_better(seeds[i], names(seeds)[i])


## =======================
## Quick sanity summary
## =======================

# summary_df <- data.frame(
#   id     = names(seeds),
#   seed   = unname(seeds),
#   k      = nchar(seeds),
#   stringsAsFactors = FALSE
# )
# Theoretical per-kmer neighbor counts (before dedup, per concrete k-mer)
# summary_df$max_per_kmer_d1 <- 3 * summary_df$k
# summary_df$max_per_kmer_d2 <- 9 * summary_df$k * (summary_df$k - 1) / 2

# If we kept per-seed results in RAM, add sizes:
# if (exists("results") && !is.null(results)) {
#   summary_df$n_unique_out <- vapply(results, length, integer(1))
# }


# View the first few lines
# head(summary_df, 10)


# degenerate_counts <- function(seed, sigma = 4L, alphabet = c("A","C","G","T")) {
#   mp <- .get_iupac_map()
#   s  <- strsplit(toupper(seed), "", fixed=TRUE)[[1]]
#   m  <- vapply(s, function(ch) length(intersect(mp[[ch]], alphabet)), integer(1))
#   E  <- prod(m)
#   N1 <- E * sum((sigma - m) / m)
#   if (length(m) >= 2) {
#     idx <- utils::combn(seq_along(m), 2L)
#     N2 <- sum((sigma - m[idx[1,]]) * (sigma - m[idx[2,]]) * E /
#                 (m[idx[1,]] * m[idx[2,]]))
#   } else N2 <- 0
#   list(k = length(s), multiplicities = m, E = as.double(E), N1 = as.double(N1), N2 = as.double(N2))
# }
# 
# sapply(as.vector(seeds), degenerate_counts)
