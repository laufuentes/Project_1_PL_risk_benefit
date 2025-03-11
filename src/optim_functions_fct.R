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

gradL <- function(psi, X, lambda, beta, centered, delta_Y, delta_Z){
  2*(psi(X) - delta_Y(X)) + lambda*sigma_beta_prime(psi,X, beta, centered)*delta_Z(X)
}

# Stochastic gradient descend 
SGD <- function(X, theta_current, lambda, beta, centered, psi, lr){
  n <- nrow(X)
  max_iter <- 1e3
  tol <- 1e-3
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  for(i in 1:max_iter){
    s <- sample(1:n, as.integer(n/3))
    x <- X[s,]
    grad <- gradL(psi, x, lambda, beta, centered, delta_Y, delta_Z)
    dstheta <-  as.matrix(2*expit(x%*%theta_current)* (1-expit(x%*%theta_current)))[,rep(1,ncol(x))] * x
    dL_dtheta <- mean(dstheta*grad) %>% as.vector()
    #if (i%%100==0){print(sum(ds_theta^2))}
    if (sum(dstheta^2) < tol) return(theta_current)
      theta_current <- as.vector(theta_current - lr * dL_dtheta)
    }
    return(theta_current)
}

# Frank-Wolve algorithm 
FW = function(X, lambda, beta, alpha, delta_Y, delta_Z){
    K <- 50
    tol <- 1e-5
    lr <-0.01
    set.seed(2025)
    theta_init <- runif(dim(X)[2], -5, 5)
    psi_terms <- list(list(weight = 1, func = function(X) { psi_theta(X, theta_init) }))
    thetas <- matrix(0,K+2, dim(X)[2])
    thetas[1,] <- theta_init
    for (k in 0:K){
      if (k%%10==0){print(k)}
        psi <- function(X) {
        # Compute the convex combination of functions and return it as a vector
        result <- sapply(psi_terms, function(p) p$weight * p$func(X))
        return(rowSums(result))  # Sum across rows to keep it as a vector
        }
        gamma <- 2/(2+k)
        theta_opt <- SGD(X, thetas[k+1,], lambda, beta, centered, psi,lr)
        thetas[k+2,] = theta_opt
        # if (k>0 & sum(abs(thetas[k+1,]-thetas[k+2,]))<tol){
        #   break
        # }
        # if (k > 0 && sum(abs(thetas[k+2,]-thetas[k+1,])) < tol) {
        #     return(list(theta = psi_terms, obj_func_vect = obj_func_vect[1:(k + 1)]))
        # }
        psi_k1 <- function(X) { psi_theta(X, theta_opt) }
    
        # Store the new psi function without deep recursion
        psi_terms <- lapply(psi_terms, function(p) {
        list(weight = (1 - gamma) * p$weight, func = p$func)})

        psi_terms[[length(psi_terms) + 1]] <- list(weight = gamma, func = psi_k1)
    }
    return(psi_terms)
}

psi <- function(X, psi_terms){
        # Compute the convex combination of functions and return it as a vector
        result <- sapply(psi_terms, function(p) p$weight * p$func(X))
        return(rowSums(result))  # Sum across rows to keep it as a vector
}
