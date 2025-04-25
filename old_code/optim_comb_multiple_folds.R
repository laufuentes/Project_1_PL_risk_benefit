args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])
if (is.na(i) || is.na(j)) stop("Need valid integer indices i and j.\n")

cat(sprintf("[START] optim_comb.R: i = %d ; j = %d\n", i, j)); flush.console()

library(parallel)
library(tidyverse)
library(grf)

# parameters
centered  <- FALSE
alpha     <- 0.05
precision <- 0.025

# Adjust these paths as needed
param_combinations<- readRDS(file.path("MC","results","data","param_combinations.rds"))

# Define or load necessary functions
source(file.path("src","algorithm_architecture.R")) 
source(file.path("src","synthetic_data.R"))

df <- read.csv(file.path("MC","results","data","estimated",paste0("df_",i,".csv")),stringsAsFactors = FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()

s <- readRDS(file.path("MC","results","data","S",paste0("s_",i,".rds")))
Jfold<- readRDS(file.path("MC","results","data","Jfold.rds"))
mu.hat.nj <- readRDS(file.path("MC","results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path("MC","results","data","Nu",paste0("nu.hat.nj_",i,".rds")))

Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
  function(X) mu.nj(1, X) - mu.nj(0, X)
})

Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
  function(X) nu.nj(1, X) - nu.nj(0, X)
})
CC_mu <-readRDS(file.path("MC","results","data","CC_mu",paste0("CC_mu_",i,".rds")))
CC_nu<-readRDS(file.path("MC","results","data","CC_nu",paste0("CC_nu_",i,".rds")))

# Run optimization for this index
for (fold in 1:Jfold){
  policy <- optimize_combination(j, X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], param_combinations, centered, alpha,precision)
  res <- parallelized_process_policy_tlearner(j, param_combinations, list(policy), X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], mu.hat.nj[[fold]], centered, alpha)
  # Save the result
  saveRDS(res, file = file.path("MC","results","estimated","individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
}
