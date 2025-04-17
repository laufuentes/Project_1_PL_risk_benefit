args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
if (is.na(i)) stop("Need a valid integer index.")

library(parallel)
library(tidyverse)

# Load or generate needed data
centered <- FALSE   # or FALSE, depending on your use case
alpha <- 0.05
precision <- 0.025

# Adjust these paths as needed
param_combinations<- readRDS(file.path("MC","results","data","param_combinations.rds"))

# Define or load necessary functions
source(file.path("src","algorithm_architecture.R")) 
source(file.path("src","synthetic_data.R"))

df_complete <- read.csv(file.path("MC","results","data","oracular",paste0("df_complete_",i,".csv")),stringsAsFactors = FALSE)

df <- read.csv(file.path("MC","results","data","estimated",paste0("df_",i,".csv")),stringsAsFactors = FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()

s <- readRDS(file.path("MC","results","data",paste0("s_",i,".rds")))
mu.hat.nj <- readRDS(file.path("MC","results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path("MC","results","data","Nu",paste0("nu.hat.nj_",i,".rds")))

Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
  function(X) mu.nj(1, X) - mu.nj(0, X)
})

Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
  function(X) nu.nj(1, X) - nu.nj(0, X)
})
CC_mu <-readRDS(file.path("MC","results","data","CC_mu",paste0("CC_mu_",i,".rds")))
CC_nu<-readRDS(file.path("MC","results","data","CC_mu",paste0("CC_nu_",i,".rds")))



# Run optimization for this index
policy <- optimize_combination(i, param_combinations,delta_mu, delta_nu,X,alpha,centered,precision)
res <- parallelized_process_policy(i, param_combinations, list(policy), X,
                      list(df_complete$y1, df_complete$y0), delta_mu, delta_nu, centered, alpha)

# Save the result
saveRDS(res, file = file.path("MC","results","oracular","individual_results",paste0("res_", i, ".rds")))

