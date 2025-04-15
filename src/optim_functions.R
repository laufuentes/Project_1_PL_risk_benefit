library(tidyverse)
library(roxygen2)

source("src/objective_functions.R")
expit <- plogis
logit <- qlogis
#' Generate Psi Function
#'
#' Constructs a convex combination function \code{psi} based on the sequence of solutions
#' obtained from the Frank-Wolfe (FW) algorithm. Each new solution \code{theta} contributes  
#' to \code{psi} in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
#'
#' @param Theta A numeric matrix of size K x d, where each row represents the k-th solution  
#'        \code{theta} obtained from the inner minimization in the FW algorithm.
#'
#' @return A function that computes \code{psi} for a given input matrix \code{x},  
#'         using a convex combination of past solutions.
#' @export
make_psi <- function(Theta) {
  # Theta: K x d
  # gamma: real
  ## ------------
  ## lazy version
  ## ------------
  psi_iterative <- function(x) {
    # x: n x d
    Theta_x <- x %*% t(Theta)
    psi_x <- 2*expit(Theta_x[,1])-1
    weights <- rep(0,K)
    weights[1] <- 1
    # Theta_x: n x K
    for (k in 1:(ncol(Theta_x)-1)){
      gamma <- 2/(2+k)
      weights[1:k] <- weights[1:k]*(1-gamma)
      weights[k+1] <- gamma
      psi_x <- (1-gamma) * psi_x + gamma * (2 * expit(Theta_x[, k+1])-1)
    }
    return(psi_x)
  }

  psi <- function(x){
    k <- nrow(Theta)
    Gamma <- matrix((2 / ((k + 1) * k)) * (1:k), ncol = 1)
    sigma_theta_x <- 2 * expit(x %*% t(Theta)) - 1 #only use the first k thetas.
    psi_x <- sigma_theta_x %*% Gamma
    return(psi_x)
  }
  return(psi)
}

#' Stochastic Gradient Descent (SGD)
#'
#' Performs stochastic gradient descent to optimize the parameters.
#'
#' @param X A numeric matrix of size n x d (input data).
#' @param theta_current A numeric matrix of size 1 x d (intialization for parameter to estimate).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to center the policy.
#' @param delta_Mu A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Nu A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#' @param psi A function that takes X as input.
#' @param verbose A logical value indicating whether to print progress.
#'
#' @return A numeric matrix of size 1 x d (optimized parameters).
#' @export
SGD <- function(X, theta_current, lambda, beta, centered, delta_Mu, delta_Nu, psi, verbose){
  n <- nrow(X)
  max_iter <- 1e3
  tol <- 1e-3
  lr <- 0.01

  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  batch_size <- as.integer(n / 3)

  LprimeX <-  gradL(psi, X, lambda, beta, centered, delta_Mu, delta_Nu)
  for(i in 1:max_iter){
    s <- sample.int(n, batch_size)
    x <- X[s,]

    theta_x <- x %*% t(theta_current)
    expit_theta_x <- expit(theta_x)
    expit_diff <- 2 * expit_theta_x * (1 - expit_theta_x)

    Lprime <-LprimeX[s]
    dL_dtheta <- t(t(x) %*% (expit_diff * Lprime))

    if (verbose && i %% 500 == 0) {
            theta_X <- X %*% t(theta_current)
            expit_theta_X_full <- expit(theta_X)
            expit_Diff <- 2 * expit_theta_X_full * (1 - expit_theta_X_full)

            Whole_Grad <- t(t(X) %*% (expit_Diff * LprimeX))

            if (mean(Whole_Grad) < tol) {
                break
            }
            value <- mean(LprimeX * (2 * expit_theta_X_full - 1))
            msg <- sprintf("SGD: iteration %i, value %f", i, value)
            message(msg)
    }

    theta_current <- theta_current - lr * dL_dtheta
    }
    return(theta_current)
}

#' Frank-Wolfe Algorithm
#'
#' Implements the Frank-Wolfe optimization algorithm to iteratively refine a convex  
#' combination function \code{psi}. At each iteration, a new solution \code{theta}  
#' is computed via stochastic gradient descent (SGD) and added to the convex combination  
#' in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
#'
#' @param X A numeric matrix (input data).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param alpha A numeric scalar (constraint tolerance).
#' @param delta_Mu A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Nu A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#' @param centered A logical value indicating whether to center the policy.
#' @param precision A numeric scalar that determines the convergence precision desired.
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the optimized parameter \code{theta},  
#'         where each row represents the k-th \code{theta} solution at iteration \code{k}.
#' @export
FW <- function(X, lambda, beta, alpha, delta_Mu, delta_Nu, centered, precision, verbose=TRUE) {
    K <- as.integer(1/precision)
    tol <- 1e-5
    d <- ncol(X)
    theta_fix <- matrix(runif(d, -5, 5), ncol=d, nrow=1)
    theta <- theta_fix
    #psi_terms <- list(list(weight = 1, func = function(X) { psi_theta(X, theta_init) }))
    
    for (k in 0:K){
      if (k==1){theta <- matrix(theta[2,], nrow=1, ncol=d)}

      psi <- make_psi(theta)
        
        if (verbose && k %% 20 == 0) {
            msg <- sprintf("FW: iteration %i, value %f", k, L(psi, X, lambda, beta, alpha, centered, delta_Mu, delta_Nu))
            message(msg)
        }
        theta_opt <- SGD(X, theta_current=theta_fix, 
                         lambda, beta, centered, delta_Mu, delta_Nu, psi, (verbose && k %% 10 == 0))
        
        theta <- rbind(theta, theta_opt)
    }
    return(theta)
}