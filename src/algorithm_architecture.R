source("src/optim_functions.R")
library(tidyverse)

#' Optimize Policy Parameters
#'
#' Computes optimized policy parameters (`theta`) for a specific combination of
#' regularization hyperparameters using the full dataset and given causal contrasts.
#'
#' @param i Integer index indicating which parameter combination to use.
#' @param param_combinations A data frame containing combinations of `lambda` and `beta` values.
#' @param delta_Y A function of estimated treatment effect contrasts (`mu`) by fold.
#' @param delta_Xi A function of estimated constraint components (`nu`) by fold.
#'
#' @return A vector of optimized policy parameters (`theta`).
#' @export
optimize_combination <- function(i, param_combinations, delta_Y, delta_Xi){
  thetas <- FW(X, 
  param_combinations$lambda[i], 
  param_combinations$beta[i], 
  alpha, delta_Y, delta_Xi, centered, precision)
  return(thetas) 
}

#' Optimize Policy Parameters Using T-Learner Estimation
#'
#' This function computes optimized policy parameters (`theta`) for a given
#' parameter combination using T-learner estimates of causal contrasts.
#'
#' @param i Integer index indicating which parameter combination to use.
#' @param param_combinations A data frame containing combinations of hyperparameters (`lambda`, `beta`) and fold IDs.
#' @param Delta_mu_nj_folds A list of estimated treatment effect contrasts (`mu`) by fold.
#' @param Delta_nu_nj_folds A list of estimated constraint components (`nu`) by fold.
#'
#' @return A vector of optimized policy parameters (`theta`).
#' @export
optimize_combination_Tlearner <- function(i, param_combinations, Delta_mu_nj_folds, Delta_nu_nj_folds){
  fold <- param_combinations$Fold[[i]]
  data <- df[s!=fold,]
  X <- data %>% select(starts_with("X.")) %>% as.matrix()
  lambda <- param_combinations$lambda[[i]]
  beta <- param_combinations$beta[[i]]
  # Step 2: Debias causal contrasts
  Delta_mu_nj <- Delta_mu_nj_folds[[fold]]
  Delta_nu_nj <- Delta_nu_nj_folds[[fold]]
  thetas <- FW(X, lambda, beta, alpha, 
  Delta_mu_nj, Delta_nu_nj, 
  centered, precision)
  return(thetas) 
}

parallelized_process_policy <- function(
     idx,
     param_combinations,
     thetas,
     covariates,
     counterfacts,
     delta_Y, delta_Xi,
     centered,
     alpha) {
     # Extract the policy for the current index
     theta <- thetas[[1]]
     psi <-make_psi(theta)
     results <- data.frame(
         lambda = param_combinations$lambda[idx],
         beta = param_combinations$beta[idx],
         optimal_x = I(list(psi(covariates))), # I() wraps the list to avoid issues with data frames
         risk = R_p(psi, covariates, delta_Y),
         constraint = S_p(
             psi,
             covariates,
             param_combinations$beta[idx],
             alpha, centered, delta_Xi
             ),
         obj = L(psi, covariates,param_combinations$lambda[idx], param_combinations$beta[idx], alpha, centered, delta_Y, delta_Xi),
         policy_value = policy_values(
           psi,
           covariates, 
           c(counterfacts[[1]], counterfacts[[2]]),
           param_combinations$beta[idx],
           centered,
           alpha
         )
     )
     colnames(results) <- c(
         "lambda",
         "beta",
         "optimal_x",
         "risk",
         "constraint",
         "obj",
         "policy_value")
     return(results) # Return the updated results for this index
 }

#' Evaluate a Policy for a Given Parameter Combination
#'
#' Evaluates a learned policy using performance metrics such as risk, constraint
#' violation, and objective value. Also computes the policy value using counterfactual outcomes.
#'
#' @param idx Index for which parameter combination and policy should be evaluated.
#' @param param_combinations Data frame of parameter combinations with columns `lambda` and `beta`.
#' @param thetas List of policy parameter vectors (output of `optimize_combination`).
#' @param covariates A matrix of covariates used to evaluate the policy.
#' @param counterfacts A list of counterfactual outcomes (e.g., `y1`, `y0`).
#' @param delta_Y A function of estimated treatment effect contrasts (`mu`) by fold.
#' @param delta_Xi A function of estimated constraint components (`nu`) by fold.
#' @param centered Logical; whether to center the features.
#' @param alpha Constraint tolerance (typically between 0 and 1).
#'
#' @return A data frame with one row containing:
#' - `lambda`, `beta`: hyperparameters
#' - `optimal_x`: list of decisions from policy
#' - `risk`: risk metric
#' - `constraint`: constraint satisfaction metric
#' - `obj`: objective value
#' - `policy_value`: estimated policy value
#' @export
process_policy <- function(
    idx,
    param_combinations,
    thetas,
    covariates,
    counterfacts,
    delta_Y, delta_Xi,
    centered,
    alpha) {
    # Extract the policy for the current index
    theta <- thetas[[idx]]
    psi <-make_psi(theta)
    results <- data.frame(
        lambda = param_combinations$lambda[idx],
        beta = param_combinations$beta[idx],
        optimal_x = I(list(psi(covariates))), # I() wraps the list to avoid issues with data frames
        risk = R_p(psi, covariates, delta_Y),
        constraint = S_p(
            psi,
            covariates,
            param_combinations$beta[idx],
            alpha, centered, delta_Xi
            ),
        obj = L(psi, covariates,param_combinations$lambda[idx], param_combinations$beta[idx], alpha, centered, delta_Y, delta_Xi),
        policy_value = policy_values(
          psi,
          covariates, 
          c(counterfacts[[1]], counterfacts[[2]]),
          param_combinations$beta[idx],
          centered,
          alpha
        )
    )
    colnames(results) <- c(
        "lambda",
        "beta",
        "optimal_x",
        "risk",
        "constraint",
        "obj",
        "policy_value")
    return(results) # Return the updated results for this index
}

#' Learn Final Policy Using Cross-Validated Contrast Estimators
#'
#' This function estimates a final policy by aggregating fold-specific contrast
#' functions (`Delta_mu` and `Delta_nu`) using cross-validation. The contrasts
#' are applied fold-wise, and the final policy parameters are estimated using
#' a policy learning algorithm (e.g., `FW`).
#'
#' @param lambda Regularization parameter for the risk component.
#' @param beta Regularization parameter for the constraint component.
#' @param X Covariate matrix (n x p) for all units.
#' @param s A vector indicating fold assignments (length n).
#' @param Delta_mu_nj_folds A list of fold-specific treatment contrast functions (e.g., T-learner estimates).
#' @param Delta_nu_nj_folds A list of fold-specific constraint contrast functions.
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
Final_policy <- function(lambda, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds){
  Delta_mu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      delta_mu <- Delta_mu_nj_folds[[fold]]
      out[s==fold] = delta_mu(X)
    }
    return(out)
  }
  Delta_nu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      delta_nu <- Delta_nu_nj_folds[[fold]]
      out[s==fold] = delta_nu(X)
    }
    return(out)
  }
  theta_final <- FW(X, lambda, beta, alpha, 
  Delta_mu_CV, Delta_nu_CV, 
  centered, precision)
}