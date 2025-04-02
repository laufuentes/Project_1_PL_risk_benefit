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
#B <- c(0.05,0.1)
#Lambda <- c(0,5,10)
B <- seq(0.05, 2, 0.05) # beta candidates 
Lambda <- seq(0,10, 1) # lambda candidates

# Estimation parameters
Jfold<- 5 # number of K-folds
technique <- "probability.forest" # technique for estimation

## Data generating process ##
exp <- data_gen(n, option)
df <- exp[[2]]
df$Xi <- ifelse(df$Z>2,1,0)

X <- df %>% select(starts_with("X")) %>% as.matrix()

param_combinations <- expand.grid(Fold=1:Jfold,lambda = Lambda, beta = B)

# Step 1: Obtain initial estimates of nuisance parameters
## Partition dataset 
s<- partition_data(Jfold, df)
## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[3]]
nu.hat.nj <- initial_nparams[[4]]

mc_cores_used <- 14

algo_combination_Tlearner <- function(comb){
    #Step 1: Choose proper parameter combination: lambda, beta and fold 
    fold <- comb$Fold
    data <- df[s!=fold,]
    X <- data %>% select(starts_with("X")) %>% as.matrix()
    lambda <- comb$lambda 
    beta <- comb$beta

    ## Select the relative nuisance parameters
    e.nj <- e.hat.nj[[fold]]
    mu.nj <- mu.hat.nj[[fold]]
    nu.nj <- nu.hat.nj[[fold]]

    # Step 2: Debias causal contrasts
    Delta_mu_nj <- function(X){mu.nj(1,X) - mu.nj(0,X)}
    Delta_nu_nj <- function(X){nu.nj(1,X) - nu.nj(0,X)}
    theta_nj <- FW(X, lambda, beta, alpha, Delta_mu_nj, Delta_nu_nj, precision)
  return(list(list(Delta_mu_nj,Delta_nu_nj), theta_nj))
}


combination_results <- mclapply(1:nrow(param_combinations), function(i) {
    if(i%%100==0){print(i)}
  algo_combination_Tlearner(param_combinations[i,])
}, mc.cores =mc_cores_used-2, mc.preschedule = TRUE)


# RowCountsDF <- data.frame(
#   CombinationID = 1:length(combination_results),  # Unique ID for each combination
#   NumRows = sapply(combination_results, function(res) nrow(res[[2]])) # Extract row count
# )

#Theta_df <- as.data.frame(do.call(rbind, lapply(combination_results, `[[`, 2)))

#save_policies_to_csv(Theta_df)
saveRDS(combination_results, file = "opt_results/estimation_T/Combination_results_TL.rds")
Causal_contrast<- lapply(combination_results, `[[`, 1)
Thetas <- lapply(combination_results, `[[`, 2)


res<-mclapply(1:nrow(param_combinations), function(i) {
    process_policy(i,param_combinations,Thetas,X,mu.hat.nj,
    Causal_contrast,
    centered,alpha)
}, mc.cores =mc_cores_used-2,mc.preschedule = TRUE)

saveRDS(res, file = "opt_results/estimation_T/res.rds")
results <- as.data.frame(do.call(rbind, res))

# Save the results to a CSV file
write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/estimation_T/", setting, ".csv")
)

result_T <- results %>%
  group_by(lambda, beta) %>%
  summarise(
    obj = mean(obj),
    risk = mean(risk),
    constraint = mean(constraint),
    policy_value = mean(policy_value),
    .groups = "drop"
  )

idx_opt_pol <- which(
  result_T$policy_value == max(result_T$policy_value[result_T$constraint <= 0])
)
idx_opt_obj <- which(result_T$obj == max(results$obj[results$constraint <= 0]))
idx_opt <- idx_opt_pol[1]

lambda_evol(
  result_T ,
  "estimation_T",option, 
  result_T$lambda[idx_opt],
  result_T$beta[idx_opt])


#### Run best policy #####
lambda <- result_T$lambda[[idx_opt]]
beta <- result_T$beta[[idx_opt]]

Final_policy <- function(lambda, beta,X, s){
  Delta_mu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      mu <- mu.hat.nj[[fold]]
      out[s==fold] = mu(1,X_fold) -mu(0,X_fold)
    }
    return(out)
  }
  Delta_nu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      nu <- nu.hat.nj[[fold]]
      out[s==fold] = nu(1,X_fold) -nu(0,X_fold)
    }
    return(out)
  }
  theta_final <- FW(X, lambda, beta, alpha, Delta_mu_CV, Delta_nu_CV, precision)
}

theta_lambda0 <- Final_policy(0, 0.05,X, s)
psi0<- make_psi(theta_lambda0)
psi_x0 <- psi0(X)

theta_opt <- Final_policy(lambda, beta,X, s)
optimal_psi <- make_psi(theta_opt)
optimal_x <- optimal_psi(X)

p0 <- gamma_plot_funct(psi_x0, 0, 0.05, df, option, centered)
ggsave("images/estimation_T/policy_lambda0.pdf",p0)
p <- gamma_plot_funct(optimal_x, result_T$lambda[[idx_opt]], result_T$beta[idx_opt], df, option, centered)
ggsave("images/estimation_T/optimal_policy.pdf",p)

#geom_points_fct(results%>%filter(fold==results$fold[[idx_opt]]), idx_opt, df, "estimation", option, centered=TRUE)

#gamma_lambda_plot(results, df, "estimation", option)


# To load them back later:
# combination_results <- readRDS("opt_results/estimation_T/Combination_results.rds")
# Thetas <- readRDS("opt_results/estimation_/Thetas.rds")
# Causal_contrast <- readRDS("opt_results/estimation/Causal_contrast.rds")
# 

res <- readRDS("opt_results/estimation_T/res.rds")
combination_results <- readRDS("opt_results/estimation_T/Combination_results_TL.rds")