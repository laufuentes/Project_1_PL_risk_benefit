set.seed(2025)

#Required libraries
library(tidyverse)
library(dplyr)

h_Y<- function(X,A, option){
  if(option=="1"){
    2*(1- X[,1]-X[,2])*A 
  }else{
    8*(1-X[,1]^2 -X[,2]*2)*A
  }
}

h_R<-function(X,A, option){
  if(option=="1"){
    (1+X[,1]-X[,2])*A
  }else{
    (X[,1]+X[,2]-0.3)*A
  }
}

data_gen <- function(n, option){
  if(option[1]=="IVF"){
    # general parameters
    u <- 10
    beta <- 3
    beta.0 <- 1
    
    theta_Y0 <- c(-0.01,0.2,-0.05, 1.5/4,2.5/4, 1.5/4,1.5/4 )
    theta_Y1<- c(0,0,0, 1, 1, 1, 1)
    
    theta_Z <- c(-0.009,0,-0.09,0,0,0,0)
    
    # Covariates 
    ## basic characteristics 
    height <- rnorm(n, mean=163, sd=7)
    age<- as.integer(runif(n, min=18,max=60))
    bmi <- rnorm(n,21.75, sd=3)
    
    ## diseases
    diabetes<- rbinom(n,1,p=ifelse(bmi>25, 0.2,4.7*1e-2))
    cancer <- rbinom(n,1, p=1e-2)
    hyperthrd <- rbinom(n,1, p=2*1e-2)
    hypothrd <- rbinom(n,1,p=ifelse(age>60,1*1e-2, 3*1e-2))
    
    X <-data.frame(height,age,bmi,diabetes, cancer, hyperthrd, hypothrd)
    
    # RCT setting
    Treatment <- rbinom(n,1,p=0.5) # Random treatment allocation
    
    # Primary outcome definition
    epsilon <- rnorm(n,mean=0, sd = 1)
    Y_0 <- as.integer(pmax(u - (as.matrix(X)%*%theta_Y0)+rnorm(n, mean=0, sd=1),0))
    Y_1 <- as.integer(pmax(Y_0 + beta*(1 - as.matrix(X)%*%theta_Y1) + rnorm(n, mean=0, sd=1),0))
    Y.0 <- (Y_0-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    Y.1 <- (Y_1-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    
    p0 <- ifelse(
      X[,3] < 20 & X[,1] < 150,
      0.5, # Generate random values for true condition
      0.055) # Generate random values for false condition
    
    p1 <- ifelse(
      X[,3] < 25 & X[,1] < 155,
      p0+ 0.35,
      p0+ 0.25)
    
    if(option[2]=="1"){
      logit_prob_Z.0 <- beta.0 + as.matrix(X)%*%theta_Z 
      logit_prob_Z.1 <- logit_prob_Z.0 + 2
      
      p0 <- 1/(1+exp(-logit_prob_Z.0))
      p1 <- 1/(1+exp(-logit_prob_Z.1))
      Z.1<- rbinom(n,1,p1)
      Z.0<- rbinom(n,1,p0)
    }
    }else{
    Treatment <- rbinom(n,1,0.5)
    X <- matrix(runif(n*10,0,1),n,10)
    epsilon_Y <- rnorm(n,0,1)
    epsilon_R <- pmin(rnorm(n,0,1),1)
    
    Y_1 <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(1,n),option[2]) + epsilon_Y
    Y_0 <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(-1,n),option[2]) + epsilon_Y
    
    Y.1 <- (Y_1-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    Y.0 <- (Y_0-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    
    p1<- 2+ X[,1] + h_R(X,rep(1,n),option[2])+epsilon_R
    p0 <- 2+ X[,1] + h_R(X,rep(-1,n),option[2])+epsilon_R
  }
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,p1=p1,p0=p0)
  df_obs<- data.frame(X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Z=ifelse(Treatment==1, p1, p0))
  return(list(df_complete, df_obs))
}

delta_Y <- function(X){
  # Access individual parameters
  if(option[1]=="IVF"){
    u <- 10
    beta <- 3
    theta_Y0 <- c(-0.01,0.2,-0.05, 1.5/4,2.5/4, 1.5/4,1.5/4 )
    theta_Y1<- c(0,0,0, 1, 1, 1, 1)
    
    Y0_X <- as.integer(
      pmax(u - (as.matrix(X)%*%theta_Y0),0)
    )
    Y1_X <- as.integer(
      pmax(
        Y0_X + beta*(1 - as.matrix(X)%*%theta_Y1),0
      )
    )
    
  }else{
    Y1_X <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(1,n),option[2])
    Y0_X <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(-1,n),option[2])
    
  }
  Y.0 <- (Y0_X-min(Y0_X,Y1_X))/(max(Y0_X,Y1_X)-min(Y0_X,Y1_X))
  Y.1 <- (Y1_X-min(Y0_X,Y1_X))/(max(Y0_X,Y1_X)-min(Y0_X,Y1_X))
  out <- Y.1 - Y.0
  return(out)
}

delta_Z <- function(X){
  theta_Z <- c(-0.009,0,-0.09,0,0,0,0)
  beta.0<- 1
  if(option[1]=="IVF"){
    if(option[2]=="1"){
      logit_prob_Z.0 <- beta.0 + as.matrix(X)%*%theta_Z 
      logit_prob_Z.1 <- logit_prob_Z.0 + 2
      
      p0 <- 1/(1+exp(-logit_prob_Z.0))
      p1 <- 1/(1+exp(-logit_prob_Z.1))
    }else{
      p0 <- ifelse(
        X[,3] < 20 & X[,1] < 150,
        0.5, # Generate random values for true condition
        0.055) # Generate random values for false condition
      
      p1 <- ifelse(
        X[,3] < 25 & X[,1] < 155,
        p0+ 0.35,
        p0+ 0.25)
    }
  }else{
    p1<- 2+ X[,1] + h_R(X,rep(1,n),option[2])
    p0 <- 2+ X[,1] + h_R(X,rep(-1,n),option[2])
  }
  out <- p1 - p0
  return(out)
}