# load_results.R
library(tidyverse)
centered <- FALSE
df <- read.csv(file.path("results","oracular","df_complete.csv"), stringsAsFactors=FALSE)

 # Directory where individual result files are saved
 res_dir <- file.path("results/oracular/individual_results")
 
 # List all result files
 res_files <- list.files(res_dir, pattern = "^res_\\d+\\.rds$", full.names = TRUE)
 
 # Load all results into a list
 all_results <- lapply(res_files, readRDS)
 
 # Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
 # This step depends on what `res` actually is in your process
 # If each `res` is a named list or a data.frame row:
 results<- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    lambda = x[[1]],
    beta = x[[2]],
    optimal_x = I(x[[3]]), # Take only the first element of optimal_x
    risk = x[[4]],
    constraint = x[[5]],
    obj = x[[6]],
    policy_value = x[[7]]
  )
}))
 
 # Save combined results
 saveRDS(results, file = "results/oracular/results.rds")
 write.csv(results%>%select(-optimal_x), file = "results/oracular/results.csv", row.names = FALSE)

 idx_opt_obj <- which(
  results$obj == max(results$obj[results$constraint <= 0])
)
idx_opt <- idx_opt_obj


source(file.path("src","plot_functions.R"))
lambda_evol(
    results, 
    "oracular", 
    results$beta[[idx_opt]]
)

geom_points_fct(results,idx_opt, df, "oracular", centered)