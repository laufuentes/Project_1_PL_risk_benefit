# run_optimization.R
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
if (is.na(i)) stop("Need a valid integer index.")

# Load required packages
library(parallel)
library(tidyverse)

# Load or generate needed data
centered <- TRUE   # or FALSE, depending on your use case
alpha <- 0.05

# Adjust these paths as needed
Lambda <- readRDS("opt_results/data/Lambda.rds")
B <- readRDS("opt_results/data/B.rds")
param_combinations <- expand.grid(lambda = Lambda, beta = B)

# Define or load necessary functions
source("src/algorithm_architecture.R") 
source("src/synthetic_data.R")
df_complete <- read.csv("opt_results/oracular/df_complete.csv",stringsAsFactors = FALSE)
X <- df_complete %>% select(starts_with("X.")) %>% as.matrix()


# Run optimization for this index
policy <- optimize_combination(i, param_combinations)
res <- process_policy(i, param_combinations, list(policy), X,
                      list(df_complete$y1, df_complete$y0), delta_Y, delta_Xi, centered, alpha)

# Save the result
saveRDS(res, file = paste0("opt_results/oracular/indiv_res/res_", i, ".rds"))
