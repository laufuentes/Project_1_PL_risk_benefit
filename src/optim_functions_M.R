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

sigma_beta <- function(psi_X, beta, centered) {
    c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered == TRUE) {
        cent <- 0.5 - c_beta * log((1 + exp(beta * 0)) / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi_X)) / (1 + exp(-beta))) + cent
    return(out)
}

# Risk function for CATE
R_p <- function(psi_X, X, delta_Y){
    out <- mean(psi_X^2 - 2 * psi_X* delta_Y(X))
    return(out)
}

# Constraint function
S_p <- function(psi_X, X, beta, alpha, centered, delta_Z){
    out <- mean(sigma_beta(psi_X, beta, centered) * delta_Z(X)) - alpha
    return(out)
}

# Objective function 
L <- function(psi_X, X,lambda, beta, alpha, centered, delta_Y, delta_Z){
    out <- R_p(psi_X, X,delta_Y) + lambda*S_p(psi_X, X, beta, alpha, centered, delta_Z)
    return(out)
}

sigma_beta_prime <- function(psi_X, beta, centered){
    c_beta <- 1 / log(
      (1 + exp(beta)) / (1 + exp(-beta))
      )
    out <- c_beta * log(
      (beta*exp(beta*psi_X))/(1+ exp(beta*psi_X))
      )
    return(out)
}


Single_step_optim <- function(X, lambda, beta, alpha, delta_Y, delta_Z){
    theta_init <- runif(dim(X)[2], -5, 5)
    psi_opt <- optim(function(psi_X){L(psi_X, X, lambda, beta, alpha, centered, delta_Y, delta_Z)}, par=2*expit(X%*%theta_init) - 1, method="L-BFGS-B", lower=-1+tol, upper=1-tol)$par
    theta_opt<- glm(logit((psi_opt +1)/2)~ X -1 )$coef %>% as.matrix()
}

# Frank-Wolve algorithm 
FW = function(X, lambda, beta, alpha, delta_Y, delta_Z){
    K <- 1e3
    tol <- 1e-5
    set.seed(2025)
    theta_init <- runif(dim(X)[2], -5, 5)
    psi_terms <- list(list(weight = 1, func = function(X) { psi_theta(X, theta_init) }))
    thetas <- matrix(0,K+2, dim(X)[2])
    thetas[1,] <- theta_init
    for (k in 0:K){
      if (k%%100==0){print(k)}
        psi <- function(X) {
        # Compute the convex combination of functions and return it as a vector
        result <- sapply(psi_terms, function(p) p$weight * p$func(X))
        return(rowSums(result))  # Sum across rows to keep it as a vector
        }
        gamma <- 2/(2+k)
        psi_opt <- optim(function(psi_X){L(psi_X, X, lambda, beta, alpha, centered, delta_Y, delta_Z)}, par=psi(X), method="L-BFGS-B", lower=-1+tol, upper=1-tol)$par
        theta_opt<- glm(logit((psi_opt +1)/2)~ X -1 )$coef %>% as.matrix()
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

