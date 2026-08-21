# =========================================================
# Expand IUPAC -> realizations (4-base), then Hamming-1
# (mutations only at non-N of the *original seed*).
# Reverse-complement canonicalisation built in; optional dedup.
# =========================================================

generate_kmers <- function(seed,
                           out_realizations = NULL,
                           out_h1          = NULL,
                           compress        = FALSE,      # pass *.gz yourself if TRUE
                           canonicalize    = TRUE,
                           dedup_realizations = FALSE,
                           dedup_h1           = FALSE,
                           log_first_h1    = 3L) {
  
  # Local maps/constants
  IUPAC_MAP <- c(
    A="A", C="C", G="G", T="T", U="T",
    R="AG", Y="CT", S="GC", W="AT", K="GT", M="AC",
    B="CGT", D="AGT", H="ACT", V="ACG",
    N="ACGT"
  )
  BASES <- c("A","C","G","T")
  
  check_iupac <- function(seq) {
    ch <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
    bad <- setdiff(unique(ch), names(IUPAC_MAP))
    if (length(bad)) stop("Found non-IUPAC characters: ", paste(bad, collapse=" "))
    ch
  }
  opts_per_pos <- function(chars) lapply(chars, function(ch) strsplit(IUPAC_MAP[[ch]], "", fixed=TRUE)[[1]])
  rev_comp <- function(s) { paste0(rev(strsplit(chartr("ACGT","TGCA", s), "", fixed = TRUE)[[1]]), collapse = "") }
  canon_one <- function(s) { if (!canonicalize) s else { rc <- rev_comp(s); if (s <= rc) s else rc } }
  
  # return character vector of H1 neighbors for one realized seq
  neighbors_h1 <- function(real_seq, mask_idx) {
    s <- strsplit(real_seq, "", fixed = TRUE)[[1]]
    L <- length(s)
    idx <- mask_idx[mask_idx >= 1L & mask_idx <= L]
    if (!length(idx)) return(character(0))
    out <- character(0)
    for (i in idx) {
      old <- s[i]
      if (is.na(old) || !nzchar(old)) next
      alts <- BASES[BASES != old]  # three bases
      if (!length(alts)) next
      # vectorised replacement for this position
      out <- c(out, vapply(alts, function(b) {
        t <- s; t[i] <- b; paste0(t, collapse = "")
      }, FUN.VALUE = character(1)))
    }
    out
  }
  
  # Streaming expansion (depth-first)
  expand_stream <- function(iupac_seq, writer) {
    chars <- check_iupac(iupac_seq); opts <- opts_per_pos(chars); k <- length(opts); buf <- character(k)
    rec <- function(pos) {
      if (pos > k) { writer(paste0(buf, collapse = "")); return(invisible(NULL)) }
      for (b in opts[[pos]]) { buf[pos] <- b; rec(pos + 1L) }
    }
    rec(1L); invisible(NULL)
  }
  
  # Prep I/O
  seedU <- toupper(seed)
  chars <- check_iupac(seedU)
  k     <- length(chars)
  mask  <- which(chars != "N")              # mutate only original non-N
  message("Mutation mask (non-N positions): ", if (length(mask)) paste(mask, collapse=",") else "<none>")
  
  conR <- conH <- NULL
  on.exit({ if (!is.null(conR)) close(conR); if (!is.null(conH)) close(conH) }, add = TRUE)
  if (!is.null(out_realizations)) conR <- if (compress) gzfile(out_realizations, "wt") else file(out_realizations, "wt")
  if (!is.null(out_h1))           conH <- if (compress) gzfile(out_h1,          "wt") else file(out_h1,          "wt")
  writeR <- if (is.null(conR)) function(x){} else function(x) writeLines(x, conR, useBytes = TRUE)
  writeH <- if (is.null(conH)) function(x){} else function(x) writeLines(x, conH, useBytes = TRUE)
  
  # Optional dedup sets
  seenR <- if (dedup_realizations) new.env(hash = TRUE, parent = emptyenv()) else NULL
  seenH <- if (dedup_h1)           new.env(hash = TRUE, parent = emptyenv()) else NULL
  seen_has <- function(env, key) !is.null(env[[key]])
  seen_add <- function(env, key) { env[[key]] <- TRUE; invisible(TRUE) }
  
  # Counters
  realN <- 0L; h1_attempt <- 0L; h1_written <- 0L; logged <- 0L
  
  # Per-realization handler (vectorised H1 core)
  handle_realization <- function(real) {
    # realizations
    r <- canon_one(real)
    if (is.null(seenR) || !seen_has(seenR, r)) {
      if (!is.null(seenR)) seen_add(seenR, r)
      writeR(r)
      realN <<- realN + 1L
    }
    # neighbors
    neigh <- neighbors_h1(real, mask)
    if (length(neigh)) {
      h1_attempt <<- h1_attempt + length(neigh)
      if (canonicalize) {
        rc <- vapply(neigh, rev_comp, FUN.VALUE = character(1))
        neigh <- ifelse(neigh <= rc, neigh, rc)
      }
      if (!is.null(seenH)) {
        keep <- vapply(neigh, function(s) { if (seen_has(seenH, s)) FALSE else { seen_add(seenH, s); TRUE } }, logical(1))
        neigh <- neigh[keep]
      }
      if (length(neigh)) {
        writeH(neigh)
        h1_written <<- h1_written + length(neigh)
        if (logged < log_first_h1) {
          nlog <- min(length(neigh), log_first_h1 - logged)
          for (i in seq_len(nlog)) message("H1 example: ", neigh[i])
          logged <<- logged + nlog
        }
      }
    }
  }
  
  # Run
  expand_stream(seedU, handle_realization)
  
  message(sprintf("Done: realizations=%s | hamming1_attempted=%s | hamming1_written=%s",
                  format(realN, big.mark=","), format(h1_attempt, big.mark=","), format(h1_written, big.mark=",")))
  invisible(list(realizations = realN, h1_attempted = h1_attempt, h1_written = h1_written, mutable_positions = mask))
}



# 1) Plain text, easy to eyeball
res1 <- generate_kmers("AGATAANN",
                       out_realizations = "/tmp/real_AGATAANN.txt",
                       out_h1          = "/tmp/h1_AGATAANN.txt",
                       compress = FALSE, canonicalize = TRUE,
                       dedup_realizations = FALSE, dedup_h1 = FALSE)
res1
# Expect: mask 1..6, realizations = 16, attempted ≈ 16*6*3 = 288, written ≈ 288 (palindromes aside)
readLines("/tmp/h1_AGATAANN.txt", n = 5)

# 2) Compressed (you pass the .gz extension yourself)
res2 <- generate_kmers("NNCCGGNNNNNNCCGGNN",
                       out_realizations = "/tmp/real_long.txt.gz",
                       out_h1          = "/tmp/h1_long.txt.gz",
                       compress = TRUE, canonicalize = TRUE,
                       dedup_realizations = FALSE, dedup_h1 = FALSE)
res2
# Then: zcat /tmp/h1_long.txt.gz | wc -l
