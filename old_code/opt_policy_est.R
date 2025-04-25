args <- commandArgs(trailingOnly = TRUE)
b <- as.integer(args[1])

beta <- 0.05
alpha <- 0.05
centered <- FALSE
precision <- 0.025

library(dplyr)
library(grf)
source(file.path("src","algorithm_architecture.R"))


l_opt <- readRDS(file.path("MC","results","estimated","lambda_opt",paste0("opt_lambda_",b,".rds")))
s <- readRDS(file.path("MC","results","data","S",paste0("s_",b,".rds")))
df <- read.csv(file.path("MC","results","data","estimated",paste0("df_",b,".csv")), stringsAsFactors=FALSE)
X <- df%>% select(starts_with("X."))%>% as.matrix()

mu.hat.nj <- readRDS(file.path("MC","results","data","Mu",paste0("mu.hat.nj_",b,".rds")))
nu.hat.nj <- readRDS(file.path("MC","results","data","Nu",paste0("nu.hat.nj_",b,".rds")))

Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
    function(X) mu.nj(1, X) - mu.nj(0, X)
})

Delta_nu_nj_folds <- lapply(nu.hat.nj, function(nu.nj) {
    function(X) nu.nj(1, X) - nu.nj(0, X)
})
theta_opt <- Final_policy(l_opt, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds, alpha, centered, precision)
saveRDS(theta_opt, file = file.path("MC","results","estimated", "theta_opt",paste0("theta_opt_",b,".rds")))

