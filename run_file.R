library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(madness)

source("src/synthetic_data.R")
source("src/estimators.R")
source("src/tool_box.R")
source("src/objective_functions.R")
source("src/optim_functions.R")

expit <- plogis
logit <- qlogis

# General parameters
n <- 1e3 # number of individuals 
alpha <- 0.1 # constraint tolerance
centered <- TRUE # centering of policy
epsilon <- 0.03 # early stopping parameter 
precision <- 0.025

##  Data generation parameters
setting <- "Other_1" # setting selected
option <- option_det(setting, "_") 

# Grid search parameters
B <- seq(0.05, 2, 0.01) # beta candidates 
Lambda <- seq(0,15, 0.1) # lambda candidates

# Estimation parameters
Jfold<- 5 # number of K-folds
technique <- "cv.glmnet" # technique for estimation

## Data generating process ##
exp <- data_gen(n, option)
df <- exp[[2]]

df$Xi <- ifelse(df$Z>2,1,0)

X <- df %>% select(starts_with("X")) %>% as.matrix()

lambda <- 1
beta <- 0.5
res <- FW(X, lambda, beta, alpha, delta_Y, delta_Z, precision)

psi<- make_psi(res)
policy_FW <- sigma_beta(psi, X, beta, centered)

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
  geom_point(aes(x=FW, y=free))+geom_abline(intercept=0,slope=1, color="red", linetype="dashed")


### ------------------------------------ ###
###  Estimation pre-required parameters  ###
### ------------------------------------ ###


c_data <- expand.grid(m=c(1,2,3),a=c(0,1))
c_data$c <- c(0,-1,-1, 1,1,1)

gdot1 <- function(X){
    rep(0, dim(as.matrix(X))[1])
}

gdot2 <-function(X){
    rep(1, dim(as.matrix(X))[1])
}

gdot3 <- function(X){
    as.matrix(rep(1,dim(as.matrix(X))[1]))
}

gdot_functs <- list(gdot1, gdot2, gdot3)

Hm<- function(a,x,m,gdot_functs,c_data, mu_nj){
    gm_dot <- gdot_functs[[m]]
    return(c_data$c[which(c_data$m==m & c_data$a==a)]*gm_dot(mu_nj(a,x)))
}

h1 <- function(a,x){
  Hm(a,x,m=1,gdot_functs,c_data, mu_nj)
}
h2 <- function(a,x){
  Hm(a,x,m=2,gdot_functs,c_data, mu_nj)
}

h_functs <- list(h1, h2)



