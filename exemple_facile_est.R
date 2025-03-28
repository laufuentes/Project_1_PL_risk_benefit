setwd("~/Project_1_PL_risk_benefit/")
library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(numDeriv)

source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/synthetic_data.R")
source("src/estimators.R")
source("src/tool_box.R")
source("src/plot_fcts_change_after.R")
source("src/algorithm_architecture.R")

######################
# General parameters #
######################
n <- 1e3 # number of individuals 
alpha <- 0.1 # constraint tolerance
centered <- TRUE # centering of policy
epsilon <- 0.03 # early stopping parameter 
precision <- 0.025

##  Data generation parameters
setting <- "Other_1" # setting selected
option <- option_det(setting, "_") 

# Grid search parameters
B <- c(0.05,0.1)
Lambda <- c(0,5,10)

# Estimation parameters
Jfold<- 3 # number of K-folds
technique <- "cv.glmnet" # technique for estimation

## Data generating process ##
exp <- data_gen(n, option)
df <- exp[[2]]
df$Xi <- ifelse(df$Z>2,1,0)

X <- df %>% select(starts_with("X")) %>% as.matrix()

#param_combinations <- expand.grid(Fold=1:Jfold,lambda = Lambda, beta = B)

# Step 1: Obtain initial estimates of nuisance parameters
## Partition dataset 
s<- partition_data(Jfold, df)
## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[3]]
nu.hat.nj <- initial_nparams[[4]]

fold <- 1
e.nj <- e.hat.nj[[fold]]
mu.nj <- mu.hat.nj[[fold]]
nu.nj <- nu.hat.nj[[fold]]

debiased_causal_contrast <- debias_procedure(df, e.nj, mu.nj, nu.nj, 0.025)

Delta_mu_nj_star <- debiased_causal_contrast[[1]]
Delta_nu_nj_star <- debiased_causal_contrast[[2]]

L(psi, X, lambda, beta, alpha, centered, delta_Y, delta_Xi)
#L(psi, X, lambda, beta, alpha, centered, Delta_mu_nj, Delta_nu_nj)
L(psi, X, lambda, beta, alpha, centered, Delta_mu_nj_star, Delta_nu_nj_star)


res_est <- FW(X, lambda, beta, alpha, Delta_mu_nj_star, Delta_nu_nj_star, precision)
psi_est<- make_psi(res_est)
policy_FW_est<- sigma_beta(psi_est, X, beta, centered)

res_oracle <- FW(X, lambda, beta, alpha, delta_Y, delta_Xi, precision)
psi_oracle<- make_psi(res_oracle)
policy_FW_oracle<- sigma_beta(psi_oracle, X, beta, centered)


tibble(
  Estimated=policy_FW_est,
  Oracle=policy_FW_oracle
) %>%
  ggplot() +
  geom_point(aes(x=Estimated, y=Oracle))+geom_abline(intercept=0,slope=1, color="red", linetype="dashed")

  ggsave("estimation_comparison.pdf")
