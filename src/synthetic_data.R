set.seed(2025)

#Required libraries
library(tidyverse)
library(dplyr)
set.seed(2025)
expit <- plogis
logit <- qlogis
#' h_Y: Treatment Effect on Y Component Function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates where each row is an observation and columns represent features.
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' h_Y(X, A)
#' @export
h_Y<- function(X,A){
  return(2*(1- X[,1]-X[,2])*A)
}

#' Synthetic Data Generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#'
#' @return A list containing two data frames: \code{df_complete} with all potential outcomes and 
#' treatment assignments, and \code{df_obs} with observed outcomes based on treatment.
#' @examples
#' data <- data_gen(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
data_gen <- function(n){
  Treatment <- stats::rbinom(n,1,0.5)
  X <- matrix(stats::runif(n*10,0,1),n,10)
  epsilon_Y <- stats::rnorm(n,0,1)

  Y.1 <- 0.5 * expit(1 - 2*X[,1] + X[,2] - X[,3] + epsilon_Y) + 
    0.5 * expit(h_Y(X,rep(1,n)))
  Y.0 <- 0.5 * expit(1 - 2*X[,1] + X[,2] - X[,3] + epsilon_Y) + 
    0.5 * expit(h_Y(X,rep(-1,n)))
    
  p1 <- expit(4*(X[,2]-1/2))
  Xi.1<- stats::rbinom(n,1,p1)
  Xi.0<- rep(0,n)
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,Xi.1=Xi.1,Xi.0=Xi.0)
  df_obs<- data.frame(X=X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs))
}

#' Conditional Average Treatment Effect Estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates.
#'
#' @return A numeric vector representing the estimated treatment effect for Y for each observation.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu(X)
#' @export
delta_mu <- function(X){
  n <- nrow(X)
  out <- 0.5*(expit(h_Y(X,rep(1,n)))-expit(h_Y(X,rep(-1,n))))
  return(out)
}


#' Conditional Average Treatment Effect Estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates.
#'
#' @return A numeric vector representing the estimated treatment effect for Xi for each observation.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu(X)
#' @export
delta_nu <- function(X){
  out <- expit(4*(X[,2]-1/2))
  return(out)
}