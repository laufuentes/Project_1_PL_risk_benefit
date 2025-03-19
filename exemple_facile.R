library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(madness)

expit <- plogis
logit <- qlogis

n <- 1e3
p <- 10
beta <- 0.05
lambda <- 2.5
alpha <- 0.1

X <- matrix(runif(n*p, -1,1),nrow=n, ncol=p) 

psi <- function(X){
    2*expit(X%*%theta)-1
}

# Probability of treating function
sigma_beta <- function(psi, X, beta, centered) {
    c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered) {
        cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi(X))) / (1 + exp(-beta))) + cent
    return(out)
}

sigma_beta_prime <- function(psi, X, beta, centered){
    if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
    c_beta <- 1 / log(
      (1 + exp(beta)) / (1 + exp(-beta))
      )
    out <- c_beta *(beta*exp(beta*psi(X)))/(1+ exp(beta*psi(X)))
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

expit_prime <- function(X){
  expit(X)*(1-expit(X))
}

# Risk function for CATE
R_p <- function(psi, X, delta_Y){
    out <- mean(psi(X)^2 - 2 * psi(X)* delta_Y(X))
    return(out)
}

# Constraint function
S_p <- function(psi, X, beta, alpha, centered, delta_Z){
    out <- mean(sigma_beta(psi,X, beta, centered) * delta_Z(X)) - alpha
    return(out)
}

# Objective function 
L <- function(psi, X,lambda, beta, alpha, centered, delta_Y, delta_Z){
    out <- R_p(psi, X,delta_Y) + lambda*S_p(psi, X, beta, alpha, centered, delta_Z)
    return(out)
}

gradL <- function(psi, X, lambda, beta, centered, delta_Y, delta_Z){
  2*(psi(X) - delta_Y(X)) + lambda*sigma_beta_prime(psi,X, beta, centered)*delta_Z(X)
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

# Frank-Wolve algorithm 
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


res <- FW(X, lambda, beta, alpha, delta_Y, delta_Z)

psi<- make_psi(res)
policy_FW <- sigma_beta(psi, X, beta, centered)

### ------------------


sigma_beta_fixed_X <- function(psi_X){
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered) {
        cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi_X)) / (1 + exp(-beta))) + cent
}

L_fixed_X <-function(psi_X){
  mean(psi_X^2 -2*psi_X*delta_Y(X) + lambda*(sigma_beta_fixed_X(psi_X)*delta_Z(X)-alpha))
}


theta <- matrix(runif(ncol(X),-1,1),nrow=1, ncol=ncol(X))
psi_opt_X <- optim(L_fixed_X, par=rep(0,n), method="L-BFGS-B", lower=-1, upper=1)$par

policy_free <- sigma_beta_fixed_X(psi_opt_X)

tibble(
  FW=policy_FW,
  free=policy_free
) %>%
  ggplot() +
  geom_point(aes(x=FW, y=free))
