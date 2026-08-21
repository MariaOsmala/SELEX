#!/usr/bin/env Rscript

suppressWarnings({
  # No external dependencies; base R only
})

library("dplyr")
library("tidyverse")

# -----------------------------
# Helpers
# -----------------------------

time_to_seconds <- function(hms) {
  if (is.na(hms) || hms == "") return(NA_real_)
  m <- regexec("^\\s*([0-9]{1,2}):([0-9]{2}):([0-9]{2})\\s*$", hms)
  parts <- regmatches(hms, m)[[1]]
  if (length(parts) != 4) return(NA_real_)
  h <- as.numeric(parts[2]); m <- as.numeric(parts[3]); s <- as.numeric(parts[4])
  return(h*3600 + m*60 + s)
}

size_to_mb <- function(s) {
  if (is.na(s) || s == "") return(NA_real_)
  m <- regexec("([0-9]+(?:\\.[0-9]+)?)\\s*(KB|MB|GB|TB)", s, ignore.case = TRUE)
  parts <- regmatches(s, m)[[1]]
  if (length(parts) != 3) return(NA_real_)
  val <- as.numeric(parts[2]); unit <- toupper(parts[3])
  factor <- switch(unit,
                   "KB" = 1/1024,
                   "MB" = 1.0,
                   "GB" = 1024.0,
                   "TB" = 1024.0*1024.0,
                   NA_real_)
  if (is.na(factor)) return(NA_real_)
  return(val * factor)
}


extract_line_value <- function(pattern, text, group = 1, ignore.case = TRUE) {
  pattern <- paste0("(?m)", pattern)            # <-- make ^/$ match line starts/ends
  m <- regexec(pattern, text, ignore.case = ignore.case, perl = TRUE)
  parts <- regmatches(text, m)[[1]]
  if (length(parts) >= group + 1) trimws(parts[group + 1]) else NA_character_
}



parse_efficiency <- function(line) {
  # e.g. "63.89% of 00:00:36 core-walltime"
  if (is.na(line) || line == "") return(list(pct=NA_real_, of=NA_character_))
  m <- regexec(":\\s*([0-9]+(?:\\.[0-9]+)?)%\\s+of\\s+(.+)$", paste0(":", line), ignore.case = TRUE, perl = TRUE)
  parts <- regmatches(paste0(":", line), m)[[1]]
  if (length(parts) != 3) return(list(pct=NA_real_, of=NA_character_))
  return(list(pct = as.numeric(parts[2]), of = trimws(parts[3])))
}

parse_job_report <- function(text) {
  # Basic fields
  job_id  <- extract_line_value("^\\s*Job ID:\\s*(.+)$", text, 1)
  array_id<- extract_line_value("^\\s*Array Job ID:\\s*(.+)$", text, 1)
  cluster <- extract_line_value("^\\s*Cluster:\\s*(.+)$", text, 1)
  usergrp <- extract_line_value("^\\s*User/Group:\\s*(.+)$", text, 1)
  state   <- toupper(ifelse(is.na(extract_line_value("^\\s*State:\\s*([A-Z_]+)", text, 1)), "", extract_line_value("^\\s*State:\\s*([A-Z_]+)", text, 1)))
  cores_s <- extract_line_value("^\\s*Cores:\\s*([0-9]+)", text, 1)
  cores   <- ifelse(is.na(cores_s), NA_integer_, suppressWarnings(as.integer(cores_s)))
  
  # Times
  cpu_util <- extract_line_value("^\\s*CPU Utilized:\\s*([0-9]{1,2}:[0-9]{2}:[0-9]{2})", text, 1)
  cpu_eff_line <- extract_line_value("^\\s*CPU Efficiency:\\s*(.+)$", text, 1)
  wall <- extract_line_value("^\\s*Job Wall-clock time:\\s*([0-9]{1,2}:[0-9]{2}:[0-9]{2})", text, 1)
  
  # Memory
  mem_used_raw <- extract_line_value("^\\s*Memory Utilized:\\s*([^\\n]+)$", text, 1)
  mem_eff_line <- extract_line_value("^\\s*Memory Efficiency:\\s*(.+)$", text, 1)
  
  # Billing
  billed_project <- extract_line_value("^\\s*Billed project:\\s*([^\\n]+)$", text, 1)
  cpu_bu <- extract_line_value("^\\s*CPU BU:\\s*([0-9.]+)", text, 1)
  mem_bu <- extract_line_value("^\\s*Mem BU:\\s*([0-9.]+)", text, 1)
  nvme_bu <- extract_line_value("^\\s*NVME BU:\\s*([0-9.]+)", text, 1)
  
  # Parse efficiencies
  ce <- parse_efficiency(cpu_eff_line)
  me <- parse_efficiency(mem_eff_line)
  
  # Normalized numerics
  cpu_util_sec <- ifelse(is.na(cpu_util), NA_real_, time_to_seconds(cpu_util))
  wall_sec <- ifelse(is.na(wall), NA_real_, time_to_seconds(wall))
  mem_used_mb <- size_to_mb(mem_used_raw)
  mem_of_mb <- size_to_mb(me$of)
  
  failed_states <- c("FAILED","CANCELLED","TIMEOUT","NODE_FAIL","OUT_OF_MEMORY","BOOT_FAIL","PREEMPTED","DEADLINE","REVOKED")
  is_running <- !is.na(state) && state == "RUNNING"
  is_failed <- !is.na(state) && state %in% failed_states
  
  # Build data.frame row
  data.frame(
    job_id = ifelse(state == "", NA, job_id),
    array_job_id = array_id,
    cluster = cluster,
    user_group = usergrp,
    state = ifelse(state == "", NA, state),
    is_running = is_running,
    is_failed = is_failed,
    cores = cores,
    cpu_utilized_hms = ifelse(is.na(cpu_util), NA, cpu_util),
    cpu_utilized_seconds = cpu_util_sec,
    cpu_efficiency_pct = ce$pct,
    cpu_efficiency_of = ifelse(is.na(ce$of), NA, ce$of),
    job_wallclock_hms = ifelse(is.na(wall), NA, wall),
    job_wallclock_seconds = wall_sec,
    memory_utilized_raw = mem_used_raw,
    memory_utilized_mb = mem_used_mb,
    memory_efficiency_pct = me$pct,
    memory_efficiency_of = ifelse(is.na(me$of), NA, me$of),
    billed_project = billed_project,
    cpu_bu = ifelse(is.na(cpu_bu), NA, as.numeric(cpu_bu)),
    mem_bu = ifelse(is.na(mem_bu), NA, as.numeric(mem_bu)),
    nvme_bu = ifelse(is.na(nvme_bu), NA, as.numeric(nvme_bu)),
    stringsAsFactors = FALSE
  )
}

parse_directory <- function(input_dir) {
  files <- paste0(input_dir, dir(input_dir))
  #fp=paste0(input_dir, "/","kmers_29346001_690.out")
  if (length(files) == 0) {
    return(data.frame())
  }
  rows <- lapply(files, function(fp){
    txt <- tryCatch(paste(readLines(fp, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
                    error = function(e) "") 
    rec <- parse_job_report(txt)
    rec$source_file <- fp
    rec
  })
  do.call(rbind, rows)
}

# -----------------------------
# Main
# -----------------------------

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 2) {
#   cat("Usage: Rscript parse_job_reports.R <input_dir> <output_csv>\n")
#   quit(status = 1)
# }

#input_dir <- args[1]
#output_csv <- args[2]




input_dir="/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/kmers_outs/"
output_csv="/projappl/project_2013895/SELEX/generate_kmers_from_seed/experiments/job_info.csv"

df <- parse_directory(input_dir)

#split array job id column

df=df%>%
  separate(array_job_id, into = c("array_job", "array_index"), sep = "_")

df <- df %>%
  mutate(across(c("job_id", "array_job", "array_index"  ),
                as.numeric))


df=df %>%
  arrange(array_job, array_index)

df <- df %>%
  mutate(rownum = row_number()) %>%
  relocate(rownum, .before = 1)  

df$array_job %>% unique()

# Write CSV
if (nrow(df) == 0) {
  warning("No matching files found or parse produced no rows.")
}
# Ensure UTF-8
con <- file(output_csv, open = "w", encoding = "UTF-8")
write.csv(df, con, row.names = FALSE, na = "")
close(con)

cat(sprintf("Wrote %d rows to %s\n", nrow(df), output_csv))
