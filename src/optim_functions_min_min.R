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


inner_element <- function(X, coefs, thetas){
    inner_element <- colSums(t(2*expit(X%*%t(thetas)) - 1) * coefs[1:nrow(coefs),])
    return(inner_element)
}

Double_min <- function(X, lambda, beta, alpha, delta_Y, delta_Z, fixed_dim, lr ){
    coef<- matrix(0, ncol=1, nrow=fixed_dim)
    theta <- matrix(0, ncol=dim(X)[2], nrow=fixed_dim)

    # Objective function 
    L(function(X)inner_element(X, coefs,thetas),X, lambda, beta, alpha, centered, delta_Y, delta_Z)

    for (i in 1:max_iter){
        coef_prev <- coef
        coef <- coef - lr*nablaL_coef(theta)
        theta <- theta - lr*nablaL_theta(coef_prev)
    }
    
}


# Frank-Wolve algorithm 
FW <- function(X, lambda, beta, alpha, delta_Y, delta_Z){
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

