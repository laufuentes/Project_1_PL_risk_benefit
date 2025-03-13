library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(roxygen2)
library(madness)

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
  psi <- function(x) {
    # x: n x d
    psi_x <- rep(0, nrow(x))
    Theta_x <- x %*% t(Theta)
    # Theta_x: n x K
    for (k in 1:ncol(Theta_x)) {
      gamma <- 2/(2+k)
      psi_x <- (1-gamma) * psi_x + gamma * (2 * expit(Theta_x[, k])-1)
    }
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
#' @param psi A function that takes X as input.
#' @param lr A numeric scalar (learning rate).
#' @param verbose A logical value indicating whether to print progress.
#'
#' @return A numeric matrix of size 1 x d (optimized parameters).
#' @export
SGD <- function(X, theta_current, lambda, beta, centered, psi, lr, verbose){
  n <- nrow(X)
  max_iter <- 1e3
  tol <- 1e-3
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  for(i in 1:max_iter){
    s <- sample(1:n, as.integer(n/3))
    x <- X[s,]
    Lprime <- gradL(psi, x, lambda, beta, centered, delta_Y, delta_Z) 
    dL_dtheta <-  t(t(x)%*%(as.matrix(2*expit(x%*%t(theta_current))* (1-expit(x%*%t(theta_current))))*Lprime) )
    if(sum(dL_dtheta^2)<tol){break}
    if (verbose) {
      if (i%%100==0){
        value <- mean(Lprime*(2*expit(x%*%t(theta_current))-1))
        msg <- sprintf("SGD: iteration %i, value %f", i, value)
        message(msg)
      } 
    } 
    #if (sum(dL_dtheta^2) < tol) {break}
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
#' @param delta_Y A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Z A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#' @param precision A numeric scalar that determines the convergence precision desired.
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the optimized parameter \code{theta},  
#'         where each row represents the k-th \code{theta} solution at iteration \code{k}.
#' @export
FW <- function(X, lambda, beta, alpha, delta_Y, delta_Z, precision, verbose=TRUE) {
    K <- as.integer(1/precision)
    tol <- 1e-5
    lr <-0.01
    theta <- matrix(runif(ncol(X), -5, 5), ncol=ncol(X), nrow=1)
    #psi_terms <- list(list(weight = 1, func = function(X) { psi_theta(X, theta_init) }))
    
    for (k in 0:K){
      if (k==1){theta <- matrix(theta[2,], nrow=1, ncol=ncol(X))}

      psi <- make_psi(theta)
      
      if (verbose) {
        if (k%%10==0){
          verbose_SGD <- TRUE
          msg <- sprintf("FW: iteration %i, value %f", k, L(psi, X,lambda, beta, alpha, centered, delta_Y, delta_Z))
          message(msg)
        }else{
          verbose_SGD <- FALSE
        } 
      } 
        theta_opt <- SGD(X,theta_current=matrix(runif(ncol(X), -5, 5), ncol=ncol(X), nrow=1) , lambda, beta, centered, psi,lr, verbose_SGD)#theta[k+1,]
        theta <- rbind(theta, theta_opt)
    }
    return(theta)
}