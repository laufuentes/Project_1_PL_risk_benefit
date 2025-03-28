setwd("~/Project_1_PL_risk_benefit/")
library(tidyverse)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)

source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/estimators.R")


# i <- 1

# for (fold in Jfold){
#     data <- df[s,]
#     X <- data %>% select(starts_with("X")) %>% as.matrix()
#     e.nj <- e.hat.nj[[fold]]
#     mu.nj <- mu.hat.nj[[fold]]
#     nu.nj <- nu.hat.nj[[fold]]
#     for (beta in B){
#         for (lambda in Lambda){
#             # Step 2: Estimation of objective function in j-th fold
#             debiased_causal_contrast <- debias_procedure(data, e.nj, mu.nj,nu.nj, lambda, beta,precision)
#             Delta_mu_nj_star <- debiased_causal_contrast[[1]]
#             Delta_nu_nj_star <- debiased_causal_contrast[[2]]

#             # Step 3: Obtain minimizer 
#             theta_nj <- FW(X, lambda, beta, alpha, Delta_mu_nj_star, Delta_nu_nj_star, precision)
#             psi <- make_psi(theta_nj)

#             # Step 4: Policy evaluation 
#        }
#     }
# }


algo_combination <- function(comb){
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
    debiased_causal_contrast <- debias_procedure(data, e.nj, mu.nj,nu.nj, lambda, beta,precision)
    Delta_mu_nj_star <- debiased_causal_contrast[[1]]
    Delta_nu_nj_star <- debiased_causal_contrast[[2]]

    # Step 3: Obtain minimizer 
    theta_nj <- FW(X, lambda, beta, alpha, Delta_mu_nj_star, Delta_nu_nj_star, precision)
  return(list(debiased_causal_contrast,theta_nj))
}

save_policies_to_csv <- function(policies, file_name = "opt_results/estimation/policies_output.csv") {
  if (is.list(policies)) {
    policies_df <- do.call(rbind, lapply(policies, as.data.frame))  # Convert list to data frame
    write.csv(policies_df, file = file_name, row.names = FALSE)
    message("Policies saved to ", file_name)
  } else {
    stop("Policies output is not a list and cannot be converted to a data frame.")
  }
}

process_policy <- function(
    idx,
    param_combinations,
    thetas,
    covariates,
    mu.hat.nj,
    debiased_causal_contrast,
    centered,
    alpha) {
    # Extract the policy for the current index
    theta <- thetas[[idx]]
    psi <-make_psi(theta)
    CC<- debiased_causal_contrast[[idx]]
    mu.nj <- mu.hat.nj[[param_combinations$Fold[[idx]]]]
    results <- data.frame(
      fold= param_combinations$Fold[idx],
        lambda = param_combinations$lambda[idx],
        beta = param_combinations$beta[idx],
        optimal_x = I(list(psi(covariates))), # I() wraps the list to avoid issues with data frames
        risk = R_p(psi, covariates, CC[[1]]),
        constraint = S_p(
            psi,
            covariates,
            param_combinations$beta[idx],
            alpha, centered, CC[[2]]
            ),
        obj = L(psi, covariates,param_combinations$lambda[idx], param_combinations$beta[idx], alpha, centered, CC[[1]], CC[[2]]),
        policy_value = policy_values(
          psi,
          covariates, 
          c(mu.nj(1,covariates), mu.nj(0,covariates)),
          param_combinations$beta[idx],
          centered,
          alpha
        )
    )
    colnames(results) <- c(
      "fold",
        "lambda",
        "beta",
        "optimal_x",
        "risk",
        "constraint",
        "obj",
        "policy_value")
    return(results) # Return the updated results for this index
}
