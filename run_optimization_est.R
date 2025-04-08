# run_optimization.R
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
if (is.na(i)) stop("Need a valid integer index.")

# Load required packages
library(parallel)
library(tidyverse)
library(grf)

# Load or generate needed data
centered <- TRUE   # or FALSE, depending on your use case
alpha <- 0.05
precision <- 0.025
technique<- "probability.forest"

# Adjust these paths as needed
Lambda <- readRDS("opt_results/data/Lambda.rds")
B <- readRDS("opt_results/data/B.rds")
Jfold <- readRDS("opt_results/data/Jfolds.rds")

param_combinations <- expand.grid(Fold=1:Jfold, lambda = Lambda, beta = B)


# Define or load necessary functions
source("src/algorithm_architecture.R") 
source("src/synthetic_data.R")
df <- read.csv("opt_results/estimation_T/df.csv",stringsAsFactors = FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()

s <- readRDS("opt_results/data/s.rds")
mu.hat.nj <- readRDS("opt_results/data/mu.hat.nj.rds")
nu.hat.nj <- readRDS("opt_results/data/nu.hat.nj.rds")

Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
  function(X) mu.nj(1, X) - mu.nj(0, X)
})

Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
  function(X) nu.nj(1, X) - nu.nj(0, X)
})
CC_mu <-readRDS("opt_results/data/CC_mu.rds")
CC_nu<-readRDS("opt_results/data/CC_nu.rds")

# Run optimization for this index
policy <- optimize_combination_Tlearner(i, param_combinations, Delta_mu_nj_folds, Delta_nu_nj_folds)
res <- parallelized_process_policy(i, param_combinations, list(policy), X,
                      c(CC_mu[[param_combinations$Fold[[i]]]],
                      CC_nu[[param_combinations$Fold[[i]]]]), 
                      Delta_mu_nj_folds[[param_combinations$Fold[[i]]]], 
                      Delta_nu_nj_folds[[param_combinations$Fold[[i]]]], 
                      centered, alpha)
                      
# Save the result
saveRDS(res, file = paste0("opt_results/estimation_T/indiv_res/res_", i, ".rds"))
