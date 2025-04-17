# Load required packages
library(parallel)
library(tidyverse)

# Load sources 
source(file.path("src","synthetic_data.R"))
source(file.path("src","estimation.R"))
source(file.path("src","cross-fitting.R"))
source(file.path("src","algorithm_architecture.R"))

# Load or generate needed data
n <- 2*2e3

exp <- data_gen(n)
df_complete <- exp[[1]]
df<- exp[[2]]
X <- df%>% select(starts_with("X."))%>% as.matrix()

centered <- TRUE   # or FALSE, depending on your use case
alpha <- 0.05
precision <- 0.025

Lambda <- seq(1,8,0.5)
beta <- 0.05
Jfold <- 5 

param_combinations <- expand.grid(lambda=Lambda, beta = beta)

folds <- CVFolds(
    N = n,
    id = NULL,
    Y = df,
    cvControl = list(V = Jfold, stratifyCV = FALSE, shuffle = TRUE)
)

folds_df <- do.call(rbind, lapply(seq_along(folds), function(v){data.frame(index=folds[[v]], fold=v)}))
s <- folds_df$fold[order(folds_df$index)]

initial_nparams <- nuissance_params(s,df)
mu.hat.nj <- initial_nparams[[2]]
nu.hat.nj <- initial_nparams[[3]]

Delta_mu_nj_folds <- vector("list", Jfold)
Delta_nu_nj_folds <- vector("list", Jfold)
CC_mu <- vector("list", Jfold)
CC_nu<- vector("list", Jfold)

# Loop over each fold to compute the contrasts
for (fold in 1:Jfold) {
  # Get the mu and nu for the current fold
  mu.nj <- mu.hat.nj[[fold]]
  nu.nj <- nu.hat.nj[[fold]]
  
  # Define the Delta_mu and Delta_nu functions for the current fold
  Delta_mu_nj_folds[[fold]] <- function(X) { mu.nj(1, X) - mu.nj(0, X) }
  Delta_nu_nj_folds[[fold]] <- function(X) { nu.nj(1, X) - nu.nj(0, X) }
  CC_mu[[fold]]<- Delta_mu_nj_folds[[fold]](X)
  CC_nu[[fold]]<-Delta_nu_nj_folds[[fold]](X)
}


theta_or <- optimize_combination(1, X, delta_mu, delta_nu, param_combinations, centered, alpha,precision)


MC_optimize_comb_est <- function(i, df,s, Delta_mu_nj_folds, Delta_nu_nj_folds, param_combinations, centered, alpha,precision){
  `%>%`<- magrittr::`%>%`
  K.folds.thetas <- list()
  for (fold in sort(unique(s))){
    data <- df[s!=fold,]
    x <- data %>% dplyr::select(dplyr::starts_with("X.")) %>% as.matrix()
    lambda <- param_combinations$lambda[[i]]
    beta <- param_combinations$beta[[i]]
    # Step 2: Debias causal contrasts
    Delta_mu_nj <- Delta_mu_nj_folds[[fold]]
    Delta_nu_nj <- Delta_nu_nj_folds[[fold]]
    thetas <- FW(x, Delta_mu_nj, Delta_nu_nj, lambda, alpha, 
                beta, centered, precision)
  K.folds.thetas<- append( K.folds.thetas, thetas)
  }
return(K.folds.thetas)
}

optimize_combination_Tlearner <- function(i, df, s, Delta_mu_nj_folds, Delta_nu_nj_folds,param_combinations, centered, alpha,precision){
  
  return(thetas) 
}