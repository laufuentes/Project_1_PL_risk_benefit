setwd("~/Project_1_PL_risk_benefit/")
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
source("src/optim_functions.R")
source("src/synthetic_data.R")
source("src/estimators.R")



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
    mu.nj <- mu.hat.nj[[fold]]
    nu.nj <- nu.hat.nj[[fold]]
    for (beta in B){

        for (lambda in L){
            # Step 2: Estimation of objective function in j-th fold
            L_nj <- Algo_debias(X, e.nj, mu.nj, nu.nj)


            # Step 3: Obtain minimizer 
            psi_nj <- algo_optim(L_nj, X, lambda, beta)

            # Step 4: Policy evaluation 

              
       }
    }
}


