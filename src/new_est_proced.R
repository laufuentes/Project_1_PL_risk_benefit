library(tidyverse)
library(grf)
source("src/optim_functions.R")
source("src/utils.R")
source("src/estimation.R")

expit <- plogis
logit <- qlogis

# alpha <- 0.1
# beta <- 0.05
# centered <- FALSE
# precision <- 0.01
# lambda <- readRDS("Simple/results/estimated/lambda_opt/opt_lambda_1.rds")

# df <- read.csv("Simple/results/data/estimated/df_1.csv",stringsAsFactors=FALSE)
# X <- df%>% select(starts_with("X."))%>% as.matrix()
# Treatment <- df$Treatment

# s <- readRDS("Simple/results/data/S/s_1.rds")
# mu <- readRDS("Simple/results/data/Mu/mu.hat.nj_1.rds")
# mu0 <- mu[[1]]
# nu <- readRDS("Simple/results/data/Nu/nu.hat.nj_1.rds")
# nu0 <- nu[[1]]

# ps.hat <- learn_propensity_score(s,X,Treatment)
# prop_score<- ps.hat[[1]]

HX <- function(a,X, prop_score){
    prop_score_X <- prop_score(1,X)
    out <- (2*a-1)/prop_score(a,X)
    return(out)
}

update_mu <- function(a,X, mu0,epsilon_1_vec,psi_vec,prop_score){
    cst <- HX(a,X, prop_score)
    out <- expit(logit(mu0(a,X))+ cst*(psi_vec%*%epsilon_1_vec))
    return(out)
}

update_nu <- function(a,X, nu0,epsilon_2_vec,sigma_psi_vec,prop_score){
    cst <- HX(a,X,prop_score)
    out <- expit(logit(nu0(a,X))+ cst*(sigma_psi_vec%*%epsilon_2_vec))
    return(out)
}


Optimization_Estimation <- function(mu0, nu0, prop_score, df, lambda, alpha, precision, beta, centered){
    tol <- 0.5
    Delta_mu <- function(X){mu0(1,X)-mu0(0,X)}
    Delta_nu <- function(X){nu0(1,X)-nu0(0,X)}

    X <- df%>% select(starts_with("X."))%>% as.matrix()
    Treatment <- df$Treatment 
    Y <- df$Y
    Xi <- df$Xi

    mu_XA <- Treatment*mu0(1,X) + (1-Treatment)*mu0(0,X)
    nu_XA <- Treatment*nu0(1,X) + (1-Treatment)*nu0(0,X)

    prop_score_X <- prop_score(1,X)
    H_cst <- HX(Treatment,X,prop_score)

    theta_init <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
    psi<- make_psi(theta_init)
    psi_X <- psi(X)

    psi_vec <- as.matrix(cbind(Inf,psi_X))
    sb <- sigma_beta(psi_X,beta, centered) 
    sigma_psi_vec <- as.matrix(sb)
    k<-0 
    while(sum((psi_vec[,ncol(psi_vec)-1] - psi_X)^2)>tol & k<5*1e2){
        new_epsilon_comb <- new_est(Y,Xi, mu_XA, nu_XA, psi_X, sb, H_cst)
        if(k==0){
            psi_vec <- psi_vec[,-1]
            epsilon_1_vec <- matrix(new_epsilon_comb[1])
            epsilon_2_vec <- matrix(new_epsilon_comb[2])
        }else{
            if(k%%10==0){
                print(psi_X)
            }
            epsilon_1_vec<- rbind(epsilon_1_vec,new_epsilon_comb[1])
            epsilon_2_vec<- rbind(epsilon_2_vec,new_epsilon_comb[2])
        }

        mu_XA <- Treatment * update_mu(1,X, mu0,epsilon_1_vec,psi_vec,prop_score) + (1-Treatment)*update_mu(0,X, mu0,epsilon_1_vec,psi_vec,prop_score)
        nu_XA <- Treatment*update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec,prop_score) + (1-Treatment)*update_nu(0,X, nu0,epsilon_2_vec,sigma_psi_vec,prop_score)

        Delta_mu <- function(X){update_mu(1,X, mu0,epsilon_1_vec,psi_vec,prop_score) - update_mu(1,X, mu0,epsilon_1_vec,psi_vec,prop_score)}
        Delta_nu <- function(X){update_nu(1,X, nu0,epsilon_2_vec,sigma_psi_vec,prop_score)-update_nu(0,X, nu0,epsilon_2_vec,sigma_psi_vec,prop_score)}

        theta <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
        psi <- make_psi(theta)
        psi_X <- psi(X)
        psi_vec <- cbind(psi_vec, psi_X)
        sb <- sigma_beta(psi(X),beta, centered)
        sigma_psi_vec <-cbind(sigma_psi_vec, sb)

        k<-k+1
    }
    return(list(epsilon_1_vec,psi_vec,epsilon_2_vec,sigma_psi_vec))
}

new_est <- function(Y, Xi, mu_prev, nu_prev, psi_new, sigma_beta_new,H_cst){
    threshold <- 0.05

    mu_XA <-  pmax(pmin(mu_prev,1-threshold),threshold)
    nu_XA <- pmax(pmin(nu_prev,1-threshold),threshold)

    mod_1 <- glm(Y~H_Y-1, offset=qlogis(mu_XA),data=data.frame(Y=Y,H_Y=H_cst*psi_new))
    mod_2 <- glm(Xi~H_Xi-1, offset=qlogis(nu_XA),data=data.frame(Xi=Xi,H_Xi=H_cst*sigma_beta_new))

    epsilon_1<- as.numeric(mod_1$coefficients)
    epsilon_2<- as.numeric(mod_2$coefficients)
    return(c(epsilon_1,epsilon_2))
}