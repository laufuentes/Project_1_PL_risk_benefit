setwd("~/Project_1_PL_risk_benefit/")
library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(fields)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(numDeriv)

source("src/objective_function.R")
source("src/optim_functions.R")
source("src/synthetic_data.R")
source("src/estimators.R")

Jfold<- 5
technique <- "cv.glmnet"

c_data <- expand.grid(m=c(1,2,3),a=c(0,1))
c_data$c <- c(0,-1,-1, 1,1,1)

gdot1 <- function(X){rep(0, dim(X)[1])}
gdot2 <-function(X){rep(1, dim(X)[1])}
gdot3 <- function(X){rep(1,dim(X)[1])}

gdot_functs <- list(gdot1, gdot2, gdot3)

Hm<- function(a,x,m,c, gdot_functs){
    gm_dot <- gdot_functs[[m]]
    return(c_data$c[which(m==m & a==a)]*gm_dot(mu_nj(a,x)))
}

# Step 1: Obtain initial estimates of nuisance parameters
## Partition dataset 
s<- partition_data(Jfold, df)
## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[3]]
nu.hat.nj <- initial_nparams[[4]]


for (fold in Jfold){
    e.nj <- e.hat.nj[[fold]]
    mu.nj <- mu1.hat.nj[[fold]]
    nu.nj <- nu1.hat.nj[[fold]]
    for (beta in B){

        for (lambda in L){
            # Step 2: Estimation of objective function in j-th fold
            Delta_mu_nj <- Algo_debias(X, e_nj, mu_nj)
            Delta_nu_nj <- Algo_debias(X, e_nj, nu_nj)

            L_nj <- function(psi, X, lambda, beta, Delta_mu_nj,Delta_nu_nj ){
                R_p(psi, X,Delta_mu_nj) + lambda*S_p(psi, X,beta, Delta_nu_nj)
            }

            # Step 3: Obtain minimizer 
            psi_nj <- algo_optim(L_nj, X, lambda, beta)

            # Step 4: Policy evaluation 

              
       }
    }
}


