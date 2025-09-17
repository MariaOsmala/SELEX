library("digest")

compare_directories <- function(dir1, dir2) {
  files1 <- list.files(dir1, full.names = FALSE)
  files2 <- list.files(dir2, full.names = FALSE)
  
  # 1. Common files
  common_files <- intersect(files1, files2)
  
  # 2. Unique to each
  only_in_dir1 <- setdiff(files1, files2)
  only_in_dir2 <- setdiff(files2, files1)
  
  cat("Files only in", dir1, ":\n")
  print(only_in_dir1)
  
  cat("\nFiles only in", dir2, ":\n")
  print(only_in_dir2)
  
  # 3. Check if common files are identical
  mismatched_files <- c()
  for (fname in common_files) {
    file1 <- file.path(dir1, fname)
    file2 <- file.path(dir2, fname)
    
    #if (!identical(readLines(file1), readLines(file2))) {
    #  mismatched_files <- c(mismatched_files, fname)
    #}
    if (digest(file1, algo = "md5", file = TRUE) != digest(file2, algo = "md5", file = TRUE)) {
      mismatched_files <- c(mismatched_files, fname)
    }
    
  }
  
  cat("\nMismatched files (different content):\n")
  print(mismatched_files)
  
  cat("\nTotal matched and identical files:", length(common_files) - length(mismatched_files), "\n")
}

# Example usage
compare_directories("/scratch/project_2013895/SELEX/kmers", "/scratch/project_2013895/SELEX/streamed_kmers")


get_file_union <- function(dir1, dir2) {
  files1 <- list.files(dir1, full.names = FALSE)
  files2 <- list.files(dir2, full.names = FALSE)
  
  # Union of filenames
  union_files <- union(files1, files2)
  
  cat("Union of files in", dir1, "and", dir2, ":\n")
  print(union_files)
  
  return(union_files)
}

# Example usage:
all_files <- get_file_union("/scratch/project_2013895/SELEX/kmers", "/scratch/project_2013895/SELEX/streamed_kmers")



