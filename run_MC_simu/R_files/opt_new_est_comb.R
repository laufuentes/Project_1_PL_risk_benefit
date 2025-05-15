args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])
fold <- as.integer(args[3])
if (is.na(i) || is.na(j) || is.na(fold)) stop("Need valid integer indices i and j.\n")

# cat(sprintf("[START] optim_comb.R: i = %d ; j = %d\n", i, j,fold)); flush.console()

library(parallel)
library(tidyverse)
library(grf)

# # parameters
centered  <- FALSE
alpha     <- 0.1
precision <- 0.025
folder <- readRDS("run_MC_simu/current_scenario.rds")
beta <- readRDS(file.path(folder,"results","data","beta.rds"))
param_combinations<- readRDS(file.path(folder,"results","data","param_combinations.rds"))

# # Define or load necessary functions
source(file.path("src","algorithm_architecture.R")) 
source(file.path(folder,"synthetic_data.R"))
source(file.path("src","new_est_proced.R"))

df <- read.csv(file.path(folder,"results","data","estimated",paste0("df_",i,".csv")),stringsAsFactors = FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()

s <- readRDS(file.path(folder,"results","data","S",paste0("s_",i,".rds")))
Jfold<- readRDS(file.path(folder,"results","data","Jfold.rds"))
mu.hat.nj <- readRDS(file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))
ps.hat <- readRDS(file.path(folder,"results","data","prop_score",paste0("ps.hat.nj_",i,".rds")))

est_res <- Optimization_Estimation(mu.hat.nj[[fold]],nu.hat.nj[[fold]],prop_score=ps.hat[[fold]],df=df, lambda=j, alpha=alpha, precision=precision, beta=beta, centered=centered)

epsilon_1_vec <- est_res[[1]]
psi_vec <- est_res[[2]]
epsilon_2_vec <- est_res[[3]]
sigma_psi_vec <- est_res[[4]]

saveRDS(est_res, file=file.path(folder,"results","new_estimated",fold,"estimation_res.rds"))

Delta_mu <- function(X){update_mu(1,X, mu0,epsilon_1_vec,psi_vec[,-ncol(psi_vec)],ps.hat[[fold]]) - update_mu(1,X, mu0,epsilon_1_vec,psi_vec[,-ncol(psi_vec)],ps.hat[[fold]])}
Delta_nu <- function(X){update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec[,-ncol(sigma_psi_vec)],ps.hat[[fold]])-update_nu(0,X, nu0,epsilon_2_vec,sigma_psi_vec[,-ncol(sigma_psi_vec)],ps.hat[[fold]])}

theta <- optimize_combination(j, X, Delta_mu, Delta_nu, param_combinations, centered, alpha,precision)
res <- parallelized_process_policy(j, param_combinations, list(theta), X, Delta_mu_nj_folds[[fold]], Delta_nu_nj_folds[[fold]], centered, alpha)
res_or <- parallelized_process_policy(j, param_combinations, list(theta), X, delta_mu, delta_nu, centered, alpha)

# # Save the result
saveRDS(res, file = file.path(folder,"results","new_estimated","individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))
saveRDS(res_or, file = file.path(folder,"results","new_estimated", "oracular","or_individual_results_per_fold",paste0("res_", i,"_",j,"_",fold,".rds")))