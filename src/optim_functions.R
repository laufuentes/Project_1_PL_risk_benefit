library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(madness)

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

# Stochastic gradient descend 
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

# Frank-Wolfe algorithm
FW <- function(X, lambda, beta, alpha, delta_Y, delta_Z, verbose=TRUE) {
    K <- 50
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
        # psi <- function(X) {
        # # Compute the convex combination of functions and return it as a vector
        # result <- sapply(psi_terms, function(p) p$weight * p$func(X))
        # return(rowSums(result))  # Sum across rows to keep it as a vector
        # }
        gamma <- 2/(2+k)
        theta_opt <- SGD(X,theta_current=matrix(runif(ncol(X), -5, 5), ncol=ncol(X), nrow=1) , lambda, beta, centered, psi,lr, verbose_SGD)#theta[k+1,]
        theta <- rbind(theta, theta_opt)
        # psi_k1 <- function(X) { psi_theta(X, theta_opt) }
    
        # Store the new psi function without deep recursion
        # psi_terms <- lapply(psi_terms, function(p) {
        # list(weight = (1 - gamma) * p$weight, func = p$func)})

        # psi_terms[[length(psi_terms) + 1]] <- list(weight = gamma, func = psi_k1)
    }
    return(theta)
}