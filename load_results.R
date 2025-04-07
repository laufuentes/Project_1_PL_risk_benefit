# load_results.R

# Directory where individual result files are saved
res_dir <- "opt_results/oracular/indiv_res"

# List all result files
res_files <- list.files(res_dir, pattern = "^res_\\d+\\.rds$", full.names = TRUE)

# Load all results into a list
all_results <- lapply(res_files, readRDS)

# Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
# This step depends on what `res` actually is in your process
# If each `res` is a named list or a data.frame row:
results <- do.call(rbind, lapply(all_results, function(x) as.data.frame(t(unlist(x)))))

# Save combined results
saveRDS(df_results, file = "opt_results/oracular/combined_results.rds")
write.csv(df_results, file = "opt_results/oracular/combined_results.csv", row.names = FALSE)