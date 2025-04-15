set.seed(2025)

#Required libraries
library(tidyverse)
library(dplyr)

h_Y<- function(X,A){
  return(2*(1- X[,1]-X[,2])*A)
}

data_gen <- function(n){
  Treatment <- rbinom(n,1,0.5)
  X <- matrix(runif(n*10,0,1),n,10)
  epsilon_Y <- rnorm(n,0,1)

  Y.1 <- 0.5 * expit(1 - 2*X[,1] + X[,2] - X[,3] + epsilon_Y) + 
    0.5 * expit(h_Y(X,rep(1,n)))
  Y.0 <- 0.5 * expit(1 - 2*X[,1] + X[,2] - X[,3] + epsilon_Y) + 
    0.5 * expit(h_Y(X,rep(-1,n)))
    
  p1 <- expit(4*(X[,2]-1/2))
  Xi.1<- rbinom(n,1,p1)
  Xi.0<- rep(0,n)
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,Xi.1=Xi.1,Xi.0=Xi.0)
  df_obs<- data.frame(X=X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs))
}

delta_mu <- function(X){
  n <- nrow(X)
  out <- 0.5*(expit(h_Y(X,rep(1,n)))-expit(h_Y(X,rep(-1,n))))
  return(out)
}

delta_nu <- function(X){
  out <- expit(4*(X[,2]-1/2))
  return(out)
}