args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])
if (is.na(i) || is.na(j)) stop("Need valid integer indices i and j.\n")

cat(sprintf("[START] optim_comb.R: i = %d ; j = %d\n", i, j)); flush.console()

library(parallel)
library(tidyverse)

# parameters
centered  <- FALSE
alpha     <- 0.05
precision <- 0.025

# Adjust these paths as needed
param_combinations<- readRDS(file.path("MC","results","data","param_combinations.rds"))

# Define or load necessary functions
source(file.path("src","algorithm_architecture.R")) 
source(file.path("src","synthetic_data.R"))

df_complete <- read.csv(file.path("MC","results","data","oracular",paste0("df_complete_",i,".csv")),stringsAsFactors = FALSE)
X <- df_complete %>% select(starts_with("X.")) %>% as.matrix()

# Run optimization for this index
theta <- optimize_combination(j, X,delta_mu, delta_nu,param_combinations,centered, alpha,precision)
res <- parallelized_process_policy(j, param_combinations, list(theta), X, delta_mu, delta_nu, centered, alpha)
# Save the result
saveRDS(theta, file = file.path("MC","results","oracular","theta_opt",paste0("res_", i,"_",j, ".rds")))
saveRDS(res, file = file.path("MC","results","oracular","individual_results",paste0("res_", i,"_",j, ".rds")))

