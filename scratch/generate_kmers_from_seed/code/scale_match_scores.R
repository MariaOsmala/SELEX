

scale_to_01 <- function(x) {
  if (length(unique(x)) == 1) {
    return(rep(1, length(x)))  # or 0, or NA — all values are the same
  }
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}