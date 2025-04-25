# load_results.R
B <- readRDS(file.path("MC","results","data","MC_iter.rds"))
Lambda <- readRDS(file.path("MC","results","data","Lambda.rds"))
L <- length(Lambda)

 # Directory where individual result files are saved
 res_dir <- file.path("MC","results", "estimated","individual_results_per_fold")

extract_j <- function(filename) {
  as.integer(sub(".*res_\\d+_(\\d+)\\.rds$", "\\1", filename))
}

for (b in 1:B){
  for (l in 1:L){
    pattern <- sprintf("^res_%d_%d_.*\\.rds$", b, l)
    res_files <- list.files(path = res_dir, pattern = pattern, full.names = TRUE)
    all_results <- lapply(res_files, readRDS)
    results<- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
        lambda = x[[1]],
        beta = x[[2]],
        risk = x[[4]],
        constraint = x[[5]],
        obj = x[[6]],
        policy_value = x[[7]])}))

    result_T <- results %>%
        summarise(
        lambda = mean(lambda),
        obj = mean(obj),
        risk = mean(risk),
        constraint = mean(constraint),
        policy_value = mean(policy_value),
        .groups = "drop")
    saveRDS(result_T, file = file.path("MC","results","estimated","individual_results",paste0("result_",b,"_",l,".rds")))      
  }
}

res_dir_ind <- file.path("MC","results", "estimated","individual_results")
extract_j <- function(filename) {
  as.integer(sub(".*result_\\d+_(\\d+)\\.rds$", "\\1", filename))
}


for (b in 1:B){
   pattern <- sprintf("result_%d_\\d+\\.rds$", b)
   res_files_ind <- list.files(path = res_dir_ind, pattern = pattern, full.names = TRUE)
    res_files_ind <- res_files_ind[order(sapply(res_files_ind, extract_j))]
    all_results <- lapply(res_files_ind, readRDS)
    results<- do.call(rbind, all_results)
    idx_opt <- which(results$lambda==min(results$lambda[results$constraint<=0]))
    saveRDS(results$lambda[[idx_opt]], file = file.path("MC","results","estimated", "lambda_opt",paste0("opt_lambda_",b,".rds")))
    write.csv(results, file = file.path("MC","results","estimated", "group_by_ind",paste0("individual_",b,".csv")))
}

for (l in 1:L){
  # List all result files
  pattern <- sprintf("^result_\\d+_%d\\.rds$", l)
  res_files_l <- list.files(path = res_dir_ind, pattern = pattern, full.names = TRUE)
  all_results <- lapply(res_files_l, readRDS)
  results<- do.call(rbind,all_results)
  saveRDS(results, file = file.path("MC","results","estimated", "group_by_lambda",paste0("lambda_",l,".rds")))
}

 # Load all results into a list
 
 # Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
 # This step depends on what `res` actually is in your process
 # If each `res` is a named list or a data.frame row:

 
 # Save combined results


# Set your folder path
folder_path <- "MC/results/estimated/individual_results_per_fold"  # Replace with your actual folder path

# List all .res files in the folder
 files <- list.files(path = folder_path, pattern = "^res_\\d+_\\d+_\\d+\\.rds$", full.names = FALSE)

# # Extract i, j, k values using regex
 matches <- regmatches(files, regexec("res_(\\d+)_(\\d+)_(\\d+)\\.rds", files))
 data <- do.call(rbind, lapply(matches, function(x) as.numeric(x[2:4])))
 colnames(data) <- c("i", "j", "k")

# # Convert to data frame
 df <- as.data.frame(data)
 library(dplyr)
# summary <- df %>%
#   group_by(i, j) %>%
#   summarise(
#     present_k = list(sort(unique(k))),
#     missing_k = list(setdiff(1:5, k)),
#     complete = length(unique(k)) == 5,
#     .groups = "drop"
#   ) %>%
#   filter(!complete)

# # Print all (i, j) combinations missing one or more k values
# print(summary, n = Inf)

# Create full expected grid of all i, j, k combinations
full <- expand.grid(i = 1:100, j = 1:17, k = 1:5)

# Merge with existing data to find missing ones
full$present <- interaction(full$i, full$j, full$k) %in% interaction(df$i, df$j, df$k)

# Find all missing rows
missing <- full[!full$present, ]

# Summarize by (i,j) which k values are missing
library(dplyr)
summary <- missing %>%
  group_by(i, j) %>%
  summarise(
    missing_k = list(k),
    .groups = "drop"
  )

# Show missing (i,j) combos and which k values are missing
print(summary, n = Inf)
