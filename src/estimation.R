library(grf)
library(tidyverse)
library(rsample)
library(glmnet)
library(roxygen2)
library(data.table)

#' Partition Data into Folds
#'
#' This function partitions a dataset into a specified number of folds for cross-validation.
#'
#' @param n.folds An integer specifying the number of folds.
#' @param df A data frame to be partitioned.
#' @return A vector of fold assignments for each observation.
#' @examples
#' partition_data(5, data.frame(x = 1:10))
#' @export
partition_data <- function(n.folds, df){
  n <- nrow(df) 
  even_split <- floor(n / n.folds) 
  if (n %% 2 != 0) { 
      s <- c(rep(1:n.folds, even_split), 1:(n - even_split * n.folds))
  }else{
    s <- c(rep(1:n.folds, even_split))
  }
      
  s <- sample(s)
  return(s)
}

#' Train Conditional Mean Outcome Models
#'
#' This function trains conditional mean outcome models for treated and control groups 
#' using `regression_forest`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param s A vector indicating the fold assignments for each observation.
#' @param X A matrix or data frame of covariates.
#' @param Treatment A binary vector (0/1) indicating treatment assignment.
#' @param Y A numeric vector or matrix of outcomes.
#' @return A list of functions predicting the expected outcome for treated units or control unis, each one trained without the j-th fold.
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' Treatment <- rbinom(100, 1, 0.5)
#' Y <- rnorm(100)
#' s <- sample(1:5, 100, replace = TRUE)
#' mu_functions <- train_cond_mean(s, X, Treatment, Y)
#' # Apply a function from the list to new data:
#' mu_functions[[1]](X[1:5, ])  # Predict treated outcomes
#' mu_functions[[2]](X[1:5, ])  # Predict control outcomes
#' }
#' @export
train_cond_mean <- function(s, X, Treatment,Y){
  n_obs <- dim(X)[1]

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  mus <- mclapply(sort(unique(s)),function(folds){
    idx <- which(s == folds) 
    Tr <- Treatment[-idx]
    XX <- X[-idx,]
    YY <- Y[-idx]
    X_nj1 <- XX[which(Tr==1),]
    X_nj0<- XX[which(Tr==0),]

    Y_nj1 <- YY[which(Tr==1)]
    Y_nj0<- YY[which(Tr==0)]
    mod_nj1<- regression_forest(X = X_nj1,Y = Y_nj1)
    mod_nj0 <- regression_forest(X = X_nj0,Y = Y_nj0)

    mu <- function(a,X){
      pred_X1 <- predict(mod_nj1,newdata = X)$predictions %>% as.vector()
      pred_X0 <- predict(mod_nj0,newdata = X)$predictions %>%as.vector()
      return(a*pred_X1 + (1-a)*pred_X0)}
    list(folds, mu)}, 
    mc.cores = detectCores(), 
    mc.preschedule = FALSE)

  mu.hat <- unlist(lapply(mus, `[[`, 2))
  return(mu.hat)
}

#' Compute Nuisance Parameters for Propensity Score and Outcome Models
#'
#' This function computes nuisance parameters such as propensity scores and 
#' conditional mean outcomes for primary (Y) and secondary (Z) for both treatment and control groups.
#'
#' @param s A vector indicating the fold assignments for each observation.
#' @param df A data frame containing the covariates, treatment assignment, primary outcome \(Y\), and the secondary outcome \( Z \).
#' @param technique A string specifying the method to use for computing propensity scores, either `"cv.glmnet"` or `"probability_forest"`.
#' @return A list containing the nuisance parameters functions and relative folds:
#'   \itemize{
#'     \item \code{s}: A vector indicating the fold assignments.
#'     \item \code{pi.hat.nj}: A list of functions predicting the propensity score for treatment assignment, trained without j-th fold.
#'     \item \code{mu.hat.nj}: A list of functions predicting the conditional mean outcome for control units, trained without j-th fold.
#'     \item \code{nu.hat.nj}: A list of functions predicting the conditional mean outcome for secondary variable \( Z \) for control units, trained without j-th fold.
#'   }
#' @examples
#' \dontrun{
#' set.seed(123)
#' df <- data.frame(X1 = rnorm(100), X2 = rnorm(100), Treatment = rbinom(100, 1, 0.5), Y = rnorm(100), Z = rbinom(100, 1, 0.5))
#' s <- sample(1:5, 100, replace = TRUE)
#' nuisance_params <- nuissance_params(s, df, "cv.glmnet")
#' }
#' @export
nuissance_params <- function(s, df, technique){
  n <- nrow(df)

  X <- df%>%select(starts_with("X."))
  Treatment <- df$Treatment
  Y <- df$Y
  Xi <- df$Xi

  #pi.hat.nj <- compute_propensity(s, X, Treatment, technique)

  mu.hat.nj <- train_cond_mean(s,X,Treatment, Y)
  nu.hat.nj<- train_cond_mean(s,X,Treatment, Xi)

  nuiss_params <- list(s=s, mu.hat.nj, nu.hat.nj)

  return(nuiss_params)
}
