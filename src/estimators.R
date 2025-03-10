library(grf)
library(rsample)
library(glmnet)
library(roxygen2)


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

  X <- df%>%select(starts_with("X"))
  Treatment <- df$Treatment
  Y <- df$Y
  Z <- df$Z

  pi.hat.nj <- compute_propensity(s, X, Treatment, technique)

  mu.hat.nj <- train_cond_mean(s,X,Treatment,Y)

  if(all(Z %in% c(0,1))){
    nu.hat.nj<- compute_propensity(s,X,Z, technique)
  }else{
    nu.hat.nj <- train_cond_mean(s, X, Treatment, Z)
  }

  nuiss_params <- list(s=s, pi.hat.nj, mu.hat.nj, nu.hat.nj)

  return(nuiss_params)
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
    X_nj1 <- X[which(Treatment[-idx]==1),]
    X_nj0<- X[which(Treatment[-idx]==0),]

    Y_nj1 <- Y[which(Treatment[-idx]==1)]
    Y_nj0<- Y[which(Treatment[-idx]==0)]
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


#' Compute Propensity Scores Using Cross-Validation
#'
#' This function estimates propensity scores using either `cv.glmnet` or `probability_forest` 
#' with cross-validation folds.
#'
#' @param s A vector indicating the fold assignments for each observation.
#' @param X A matrix or data frame of covariates.
#' @param Treatment A binary vector (0/1) indicating treatment assignment.
#' @param technique A string specifying the method to use, either `"cv.glmnet"` or `"probability_forest"`.
#' @return A list of functions, each computing propensity scores for `a,X`, trained without the j-th fold.
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' Treatment <- rbinom(100, 1, 0.5)
#' s <- sample(1:5, 100, replace = TRUE)
#' compute_propensity(s, X, Treatment, "cv.glmnet")
#' }
#' @export
compute_propensity <- function(s, X, Treatment, technique){
  n_obs <- dim(X)[1]
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if(technique=="cv.glmnet"){
    ps <- mclapply(sort(unique(s)), function(folds){
      idx <- which(s == folds)    
      X_nj <- X[-idx,]
      T_nj <- Treatment[-idx]
      mod_nj <- cv.glmnet(x = X_nj, y = T_nj, family = "binomial")
      # Predict on left-out fold
      preds <- function(a,X){
        pred_X <- predict(mod_nj,newx = X, type = "response", s = "lambda.1se")
        return(a*pred_X + (1-a)*(1-pred_X))}  
      list(folds, preds)}, mc.cores = detectCores(), mc.preschedule = FALSE)

  }else{
    ps <- mclapply(unique(s), function(folds){
      idx <- which(s == folds)
      X_nj <- X[-idx,]
      T_nj <- as.factor(Treatment[-idx])
      mod_nj <-probability_forest(X_nj,T_nj)
      preds <- function(a,X){
         pred_X<- predict(mod_nj ,newdata = X,type = "response")$pred[, 2]
        return(a*pred_X + (1-a)*(1-pred_X))}
        list(idx, preds)
      }, mc.cores = detectCores(), mc.preschedule = FALSE)
      }
  pi.hat_nj <- unlist(lapply(ps, `[[`, 2))
  return(pi.hat_nj)
}


# dr_learner <- function(df, primary){
#   n <- nrow(df)
#   X <- df%>%select(starts_with("X"))
#   Treatment <- df$Treatment
#   Y <- df$Y

#   tau_DR <- rep(0, n)
#   s<- df$s
#   pi.hat<- df$pi.hat
#   if(primary==TRUE){
#     mu1.hat<- df$mu1.hat
#     mu0.hat<- df$mu0.hat
#   }else{
#     mu1.hat<- df$nu1.hat
#     mu0.hat<- df$nu0.hat
#   }
  

#   if (!is.matrix(X)) {
#     X <- as.matrix(X)
#   }

#   if (!is.matrix(Y)) {
#     Y <- as.matrix(Y)
#   }

#   pseudo <- (
#     (
#       Treatment - pi.hat
#       ) / (
#         pi.hat * (1 - pi.hat)
#         )) * 
#         (
#           Y - Treatment * mu1.hat - (1 - Treatment) * mu0.hat
#           ) + mu1.hat - mu0.hat

#   dr_lr <- mclapply(unique(s), function(fold){
#     idx <- which(s==fold)

#   tau.hat <- predict(
#     cv.glmnet(X[-idx, ], pseudo[-idx]), 
#     newx = X[idx,], 
#     s = "lambda.min")
    
#     list(idx, tau.hat)
#   }, 
#   mc.cores = detectCores(), 
#   mc.preschedule = FALSE)

#   tau_DR<- tau_DR[unlist(lapply(dr_lr, `[[`, 1))] <- unlist(lapply(dr_lr, `[[`, 2))
#   return(tau_DR)
# }


debias_procedure <- function(X, outcome, h_functs, e_nj, f_nj){
  varphi_k <- function(x) x

  # 1- Construct varphi_hat
  varphi_hat <- function(a,x){
    h_functs[[1]]*varphi_k(x)
    h_functs[[2]]*varphi_k(x)
  }

  # 2- Choose link function and obtain regresion coefficients
  if(all(outcome %in% c(0,1))){
    g <- function(x){logit(x)}
    g_bar <- function(x){expit(x)}

    mod <- glm(outcome ~ varphi_hat(a,x), offset=g(f_nj(a,x)), weight=(1/e_nj(a,x)))
  }
  else{
    g<- function(x){x}
    g_bar <- function(x){1/x}
  }
  # Return debiased out-of fold outcome regression estimator
  f_nj.star <- function(a, x){g_bar(g(f_nj(a,x)) + t(varphi_hat(a,x))%*%beta_nj)}
}