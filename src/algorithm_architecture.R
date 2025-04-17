source(file.path("src","optim_functions.R"))
library(tidyverse)
#' Optimize Policy Parameters
#'
#' Computes optimized policy parameters (`theta`) for a specific combination of
#' regularization hyperparameters using the full dataset and given causal contrasts.
#'
#' @param i Integer index indicating which parameter combination to use.
#' @param X Covariate matrix (n x p) for all units.
#' @param delta_Mu A function of estimated treatment effect contrasts (`mu`) by fold.
#' @param delta_Nu A function of estimated constraint components (`nu`) by fold.
#' @param param_combinations A data frame containing combinations of `lambda` and `beta` values.
#' @param centered Logical; whether to center the sigma_beta
#' @param alpha Constraint tolerance (typically between 0 and 1).
#' @param precision A numeric scalar that determines the convergence precision desired.
#'
#' @return A vector of optimized policy paameters (`theta`).
#' @export
optimize_combination <- function(i, X, delta_Mu, delta_Nu, param_combinations, centered, alpha,precision){
  thetas <- FW(X, delta_Mu, delta_Nu, 
  param_combinations$lambda[i], alpha,
  param_combinations$beta[i], centered, precision)
  return(thetas) 
}

#' Optimize Policy Parameters Using T-Learner Estimation
#'
#' This function computes optimized policy parameters (`theta`) for a given
#' parameter combination using T-learner estimates of causal contrasts.
#'
#' @param i Integer index indicating which parameter combination to use.
#' @param df A data frame containing the covariates, treatment assignment, primary outcome \(Y\), and the secondary outcome \( Xi \).
#' @param s A vector indicating the fold assignments for each observation.
#' @param Delta_mu_nj_folds A list of estimated treatment effect contrasts (`mu`) by fold.
#' @param Delta_nu_nj_folds A list of estimated constraint components (`nu`) by fold.
#' @param param_combinations A data frame containing combinations of hyperparameters (`lambda`, `beta`) and fold IDs.
#' @param centered Logical; whether to center the sigma_beta.
#' @param alpha Constraint tolerance (typically between 0 and 1).
#' @param precision A numeric scalar that determines the convergence precision desired.
#'
#' @return A vector of optimized policy parameters (`theta`).
#' @export
optimize_combination_Tlearner <- function(i, df, s, Delta_mu_nj_folds, Delta_nu_nj_folds,param_combinations, centered, alpha,precision){
  `%>%`<- magrittr::`%>%`
  fold <- param_combinations$Fold[[i]]
  data <- df[s!=fold,]
  X <- data %>% dplyr::select(dplyr::starts_with("X.")) %>% as.matrix()
  lambda <- param_combinations$lambda[[i]]
  beta <- param_combinations$beta[[i]]
  # Step 2: Debias causal contrasts
  Delta_mu_nj <- Delta_mu_nj_folds[[fold]]
  Delta_nu_nj <- Delta_nu_nj_folds[[fold]]
  thetas <- FW(X, Delta_mu_nj, Delta_nu_nj, lambda, alpha, 
               beta, centered, precision)
  return(thetas) 
}

#' Evaluate one result policy for a Given Parameter Combination
#'
#' Evaluates a learned policy using performance metrics such as risk, constraint
#' violation, and objective value. Also computes the policy value using counterfactual outcomes.
#'
#' @param idx Index for which parameter combination should be evaluated.
#' @param param_combinations Data frame of parameter combinations with columns `lambda` and `beta`.
#' @param thetas List of one policy parameter vectors (output of `optimize_combination`).
#' @param X A matrix of covariates used to evaluate the policy.
#' @param counterfacts A list of counterfactual outcomes (in order: `y1`, `y0`).
#' @param delta_Mu A function of estimated treatment effect contrasts (`mu`) by fold.
#' @param delta_Nu A function of estimated constraint components (`nu`) by fold.
#' @param centered Logical; whether to center the sigma_beta
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
parallelized_process_policy <- function(
     idx,
     param_combinations,
     thetas,
     X,
     counterfacts,
     delta_Mu, delta_Nu,
     centered,
     alpha) {
     # Extract the policy for the current index
     theta <- thetas[[1]]
     psi <-make_psi(theta)
     results <- data.frame(
         lambda = param_combinations$lambda[idx],
         beta = param_combinations$beta[idx],
         optimal_x = I(list(psi(X))), # I() wraps the list to avoid issues with data frames
         risk = R_p(psi, X, delta_Mu),
         constraint = S_p(
             psi,
             X,
             param_combinations$beta[idx],
             alpha, centered, delta_Nu
             ),
         obj = Lagrangian_p(psi, X,param_combinations$lambda[idx], param_combinations$beta[idx], alpha, centered, delta_Mu, delta_Nu),
         policy_value = V_p(
           psi,
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
#' @param X A matrix of covariates used to evaluate the policy.
#' @param counterfacts A list of counterfactual outcomes (e.g., `y1`, `y0`).
#' @param delta_Mu A function of estimated treatment effect contrasts (`mu`) by fold.
#' @param delta_Nu A function of estimated constraint components (`nu`) by fold.
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
    X,
    counterfacts,
    delta_Mu, delta_Nu,
    centered,
    alpha) {
    # Extract the policy for the current index
    theta <- thetas[[idx]]
    psi <-make_psi(theta)
    results <- data.frame(
        lambda = param_combinations$lambda[idx],
        beta = param_combinations$beta[idx],
        optimal_x = I(list(psi(X))), # I() wraps the list to avoid issues with data frames
        risk = R_p(psi, X, delta_Mu),
        constraint = S_p(
            psi,
            X,
            param_combinations$beta[idx],
            alpha, centered, delta_Nu
            ),
        obj = Lagrangian_p(psi, X, delta_Mu, delta_Nu, param_combinations$lambda[idx], alpha, param_combinations$beta[idx], centered),
        policy_value = V_p(
          psi,
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
#' @param alpha Constraint tolerance (typically between 0 and 1).
#' @param centered Logical; whether to center the features.
#' @param precision A numeric scalar that determines the convergence precision desired.
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
Final_policy <- function(lambda, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds, alpha, centered, precision){
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
  theta_final <- FW(X, Delta_mu_CV, Delta_nu_CV, lambda, alpha, 
                    beta, centered, precision)
  return(theta_final)
}