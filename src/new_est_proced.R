library(tidyverse)
library(grf)
source("src/optim_functions.R")
source("src/utils.R")
source("src/estimation.R")

expit <- plogis
logit <- qlogis

alpha <- 0.1
beta <- 0.05
centered <- FALSE
precision <- 0.01
lambda <- readRDS("Simple/results/estimated/lambda_opt/opt_lambda_1.rds")

df <- read.csv("Simple/results/data/estimated/df_1.csv",stringsAsFactors=FALSE)
X <- df%>% select(starts_with("X."))%>% as.matrix()
Treatment <- df$Treatment

s <- readRDS("Simple/results/data/S/s_1.rds")
mu <- readRDS("Simple/results/data/Mu/mu.hat.nj_1.rds")
mu0 <- mu[[1]]
nu <- readRDS("Simple/results/data/Nu/nu.hat.nj_1.rds")
nu0 <- nu[[1]]

ps.hat <- learn_propensity_score(s,X,Treatment)
prop_score<- ps.hat[[1]]

HX <- function(a,X){
    prop_score_X <- prop_score(1,X)
    out <- (a/prop_score_X)- ((1-a)/(1-prop_score_X))
    return(out)
}

update_mu <- function(a,X, mu0,epsilon_1_vec,psi_vec){
    cst <- HX(a,X)
    out <- expit(logit(mu0(a,X))+ cst*(psi_vec%*%epsilon_1_vec))
    return(out)
}

update_nu <- function(a,X, nu0,epsilon_2_vec,sigma_psi_vec){
    cst <- HX(a,X)
    out <- expit(logit(nu0(a,X))+ cst*(sigma_psi_vec%*%epsilon_2_vec))
    return(out)
}


Optimization_Estimation <- function(mu0, nu0, prop_score,df){
    tol <- 1e-2
    Delta_mu <- function(X){mu0(1,X)-mu0(0,X)}
    Delta_nu <- function(X){nu0(1,X)-nu0(0,X)}

    X <- df%>% select(starts_with("X."))%>% as.matrix()
    prop_score_X <- prop_score(1,X)

    theta_init <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
    psi<- make_psi(theta_init)
    psi_vec <- psi(X)
    sigma_psi_vec <- sigma_beta(psi_vec,beta, centered)
    epsilon_vec <- matrix(c(0,0),nrow=1,ncol=2)

    while(sum((psi_vec[,ncol(psi_vec)] - psi_X)^2)>tol){
        epsilon_vec <- new_est(df, mu0,mu0,prop_score_X, psi_vec,sigma_psi_vec, epsilon_vec)
        epsilon_1_vec<- epsilon_vec[,1]
        epsilon_2_vec<- epsilon_vec[,2]
        Delta_mu <- function(X){update_mu(1,X, mu0,epsilon_1_vec,psi_vec) - update_mu(1,X, mu0,epsilon_1_vec,psi_vec)}
        Delta_nu <- function(X){update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec)-update_nu(0,X, nu0,epsilon_2_vec,sigma_psi_vec)}

        theta <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
        psi <- make_psi(theta)
        psi_X <- psi(X)
        psi_vec <- cbind(psi_vec, psi_X)
        sigma_psi_vec <-cbind(sigma_psi_vec, sigma_beta(psi(X),beta, centered))
    }
}

new_est <- function(df, mu0,nu0,prop_score_X, psi_vec,sigma_psi_vec, epsilon_vec){
    threshold <- 0.05
    X <- df %>% select(starts_with("X."))%>%as.matrix()
    Treatment <- df$Treatment
    Y <- df$Y
    Xi <- df$Xi
    H_cst <- HX(Treatment,X)

    epsilon_1_vec<- epsilon_vec[,1]
    epsilon_2_vec<- epsilon_vec[,2]
    mu_XA <-  pmax(pmin(Treatment * update_mu(1,X, mu0,epsilon_1_vec,psi_vec) + (1-Treatment)*update_mu(0,X, mu0,epsilon_1_vec,psi_vec),1-threshold),threshold)
    nu_XA <- pmax(pmin(Treatment * update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec) + (1-Treatment)*update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec),1-threshold),threshold)

    mod_1 <- glm(Y~H_Y-1, offset=qlogis(mu_XA),data=data.frame(Y=Y,H_Y=H_cst*psi_vec[,ncol(psi_vec)]))
    mod_2 <- glm(Xi~H_Xi-1, offset=qlogis(nu_XA),data=data.frame(Xi=Xi,H_Xi=H_cst*sigma_psi_vec[,ncol(sigma_psi_vec)]))

    epsilon_1<- as.numeric(mod_1$coefficients)
    epsilon_2<- as.numeric(mod_2$coefficients)
    epsilon_vec <- rbind(epsilon_vec,c(epsilon_1,epsilon_2))
    return(epsilon_vec)
}

