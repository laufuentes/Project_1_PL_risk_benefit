library(tidyverse)

process_results <- function(res_files) {
  # Read the files and extract required information
  all_results <- lapply(res_files, readRDS)
  
  # Combine the results into a single data frame
  do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      lambda = x[[1]],
      beta = x[[2]],
      risk = x[[4]],
      constraint = x[[5]],
      obj = x[[6]],
      policy_value = x[[7]]
    )
  }))
}

summarize_results <- function(results) {
  # Summarize the results (calculate means)
  results %>%
    summarise(
      lambda = mean(lambda),
      obj = mean(obj),
      risk = mean(risk),
      constraint = mean(constraint),
      policy_value = mean(policy_value),
      .groups = "drop"
    )
}

extract_i <- function(filename) {
  as.integer(sub(".*res_(\\d+)_\\d+\\.rds$", "\\1", filename))
}

extract_j <- function(filename) {
  as.integer(sub(".*res_\\d+_(\\d+)\\.rds$", "\\1", filename))
}

extract_j_csv <- function(filename) {
  as.integer(sub(".*result_\\d+_(\\d+)\\.csv$", "\\1", filename))
}