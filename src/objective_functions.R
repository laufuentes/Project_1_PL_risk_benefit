library(tidyverse)
library(parallel)

# Probability of treating function
sigma_beta <- function(psi, X, beta, centered) {
    c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered == TRUE) {
        cent <- 0.5 - c_beta * log((1 + exp(beta * 0)) / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi(X))) / (1 + exp(-beta))) + cent
    return(out)
}

# Derivative of sigma_beta with respect to psi
sigma_beta_prime <- function(psi, X, beta, centered){
    if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
    c_beta <- 1 / log(
      (1 + exp(beta)) / (1 + exp(-beta))
      )
    out <- c_beta * log(
      (beta*exp(beta*psi(X)))/(1+ exp(beta*psi(X)))
      )
    return(out)
}

# Extrema point 
psi_theta <- function(X, theta){
    if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  out <- 2*expit(X%*%theta) -1
  return(out)
}

# Risk function for CATE
R_p <- function(psi, X, delta_Y){
    out <- mean(psi(X)^2 - 2 * psi(X)* delta_Y(X))
    return(out)
}

# Constraint function
S_p <- function(psi, X, beta, alpha, centered, delta_Z){
    out <- mean(sigma_beta(psi, beta, centered) * delta_Z(X)) - alpha
    return(out)
}

# Objective function 
L <- function(psi, X,lambda, beta, alpha, delta_Y, delta_Z){
    out <- R_p(psi, X,delta_Y) + lambda*S_p(psi, X, beta, alpha, centered, delta_Z)
    return(out)
}

policy_values <- function(psi, X,counterfactual_outcomes,beta,centered,alpha){
  sigma_psi <- sigma_beta(psi, X, beta, centered)
  policy <- rbinom(nrow(X), 1, sigma_psi)
  y1 <- counterfactual_outcomes[1]
  y0 <- counterfactual_outcomes[2]
  out <- mean(policy * y1 + (1 - policy) * y0)
  return(out)
}