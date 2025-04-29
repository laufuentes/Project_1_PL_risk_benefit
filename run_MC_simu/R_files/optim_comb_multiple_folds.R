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
alpha     <- 0.1
precision <- 0.025
folder <- readRDS("run_MC_simu/current_scenario.rds")

# Adjust these paths as needed
param_combinations<- readRDS(file.path(folder,"results","data","param_combinations.rds"))

# Define or load necessary functions
source(file.path("src","algorithm_architecture.R")) 
source(file.path(folder,"synthetic_data.R"))

df <- read.csv(file.path(folder,"results","data","estimated",paste0("df_",i,".csv")),stringsAsFactors = FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()

s <- readRDS(file.path(folder,"results","data","S",paste0("s_",i,".rds")))
Jfold<- readRDS(file.path(folder,"results","data","Jfold.rds"))
mu.hat.nj <- readRDS(file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))

Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
  function(X) mu.nj(1, X) - mu.nj(0, X)
})

Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
  function(X) nu.nj(1, X) - nu.nj(0, X)
})

# Run optimization for this index
for (fold in 1:Jfold){
  theta <- optimize_combination(j, X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], param_combinations, centered, alpha,precision)
  res <- parallelized_process_policy(j, param_combinations, list(theta), X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], centered, alpha)
  res_or <- parallelized_process_policy(j, param_combinations, list(theta), X, delta_mu, delta_nu, centered, alpha)
  # Save the result
  saveRDS(res, file = file.path(folder,"results","estimated","individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
  saveRDS(res_or, file = file.path(folder,"results","estimated", "oracular","or_individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
}

# B <- readRDS(file.path(folder,"results","data","MC_iter.rds"))
# Lambda <- readRDS(file.path(folder,"results","data","Lambda.rds"))
# L <- length(Lambda)
# Jfold<- readRDS(file.path(folder,"results","data","Jfold.rds"))


# for (i in 1:B){
#   for (j in 1:L){
#     df <- read.csv(file.path(folder,"results","data","estimated",paste0("df_",i,".csv")),stringsAsFactors = FALSE)
#     X <- df %>% select(starts_with("X.")) %>% as.matrix()
    
#     s <- readRDS(file.path(folder,"results","data","S",paste0("s_",i,".rds")))
#     mu.hat.nj <- readRDS(file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
#     nu.hat.nj <- readRDS(file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))

#     Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
#       function(X) mu.nj(1, X) - mu.nj(0, X)
#     })

#     Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
#       function(X) nu.nj(1, X) - nu.nj(0, X)
#     })

#     # Run optimization for this index
#     for (fold in 1:Jfold){
#       theta <- optimize_combination(j, X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], param_combinations, centered, alpha,precision)
#       res <- parallelized_process_policy(j, param_combinations, list(theta), X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], centered, alpha)
#       res_or <- parallelized_process_policy(j, param_combinations, list(theta), X, delta_mu, delta_nu, centered, alpha)
#       # Save the result
#       saveRDS(res, file = file.path(folder,"results","estimated","individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
#       saveRDS(res_or, file = file.path(folder,"results","estimated","or_individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
#     }
#   }
# }