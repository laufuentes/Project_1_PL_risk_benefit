library(grf)
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
  Z <- df$Xi

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

# J_funct_appprox <- function(a,e_xa, psi_x, lambda, sb) {  
#   # Compute the inner term efficiently
#   inner <- ((2*a-1) / e_xa) * (psi_x + lambda * sb)^2
  
#   return(mean(inner))  # Mean operation remains the same
# }

# gradj <- function(a,e_xa, psi_x, lambda, sb, sb_prime) {  
#   # Compute the gradient efficiently
#   term <- 2*(psi_x + lambda*sb)*(1+lambda*sb_prime)
  
#   return(( (2*a-1)/ e_xa) * (term))
# }

J_funct_appprox <- function(e_xa, psi_x, lambda, rho11, rho22, sb){  
  inner <- -2*(rho11*psi_x) +psi_x^2 + lambda*(-2*rho22*sb +sb^2)
  return(mean((1/e_xa)*inner))  # Mean operation remains the same
}

gradj <- function(e_xa,rho1, rho2){
  out <- rho1+rho2
    return( (1/e_xa) *out)
}



#' Stochastic Gradient Descent (SGD)
#'
#' Performs stochastic gradient descent to optimize the parameters.
#'
#' @param df A data.frame n x d (input data).
#' @param theta_current A numeric matrix of size 1 x d (intialization for parameter to estimate).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to center the policy.
#' @param psi A treatment rule function that takes X as input.
#' @param e_n The propensity score function, taking A and X as input.
#' @param mu_n The conditional expectation of Y function, taking A and X as input.
#' @param nu_n The conditional expectation of Xi function, taking A and X as input.
#' @param verbose A logical value indicating whether to print progress.
#'
#' @return A numeric matrix of size 1 x d (optimized parameters).
#' @export
SGD_estimation <- function(df, theta_current, lambda, beta, centered, psi, e_n, mu_n, nu_n, verbose){
  n <- nrow(df)
  x_cols <- grep("^X", names(df))
  setDT(df)
  max_iter <- 5*1e3
  tol <- 1e-3
  lr <- 0.01

  A <- df$Treatment
  X <- as.matrix(df[, ..x_cols])
  e_XA <- e_n(A, X)
  mu_XA <- mu_n(A, X)
  nu_XA <- nu_n(A, X)
  psi_X <- psi(X)
  sbX <- sigma_beta(psi,X, beta, centered)
  sb_primeX <- sigma_beta_prime(psi,X, beta)

  batch_size <- as.integer(n / 3)

  for(i in 1:max_iter){
    s <- sample.int(n, batch_size, replace = FALSE)
    data <- df[s]
    x <- as.matrix(data[, ..x_cols])
    a <- data[["Treatment"]]
    y <- data[["Y"]]
    xi <- data[["Xi"]]

    # Precompute values that don't depend on theta
    psi_x <- psi_X[s]
    mu_xa <- mu_XA[s]
    nu_xa <- nu_XA[s]
    e_xa<- e_XA[s]
    sb <- sbX[s]
    sb_prime <- sb_primeX[s]
    a_scaled <- 2 * a - 1

    # Calculate e_xa, mu_xa, and nu_xa

    rho1 <- psi_x - (y - mu_xa) * a_scaled
    rho2 <- lambda * sb_prime * (sb - (xi - nu_xa) * a_scaled)
    Jprime <- (rho1 + rho2) / e_xa

    theta_x <- x %*% t(theta_current)
    expit_theta_x <- expit(theta_x)
    expit_diff <- 2 * expit_theta_x * (1 - expit_theta_x)

    dJ_dtheta <- t(x) %*% (expit_diff * Jprime) / batch_size
    dJ_dtheta <- t(t(x) %*% (expit_diff * Jprime))

    if (verbose && i %% 1000 == 0) {
            if (mean(dJ_dtheta) < tol) {
                break
            }
            value <- mean(Jprime * (2 * expit_theta_x - 1))
            msg <- sprintf("SGD-Est: iteration %i, value %f", i, value)
            message(msg)
    }

    theta_current <- theta_current - lr * dJ_dtheta
    }
    return(theta_current)
}

#' Frank-Wolfe Algorithm for Estimation
#'
#' Implements the Frank-Wolfe optimization algorithm to iteratively refine a convex  
#' combination function \code{psi}. At each iteration, a new solution \code{theta}  
#' is computed via stochastic gradient descent (SGD) and added to the convex combination  
#' in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
#'
#' @param df A data frame n x d (input data).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to center the policy.
#' @param e_n The propensity score function, taking A and X as input.
#' @param mu_n The conditional expectation of Y function, taking A and X as input.
#' @param nu_n The conditional expectation of Xi function, taking A and X as input.
#' @param precision A numeric scalar that determines the convergence precision desired.
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the optimized parameter \code{theta},  
#'         where each row represents the k-th \code{theta} solution at iteration \code{k}.
#' @export
FW_estimation <- function(df, lambda,  beta, centered, e_n, mu_n, nu_n,precision, verbose=TRUE) {
    K <- as.integer(1/precision)
    tol <- 1e-5

    Y <- matrix(df$Y)
    X <- df %>% select(starts_with("X")) %>% as.matrix()
    A <- matrix(df$Treatment)
    Xi<- matrix(ifelse(df$Z==1,1,0))
    d <- ncol(X)

    theta_fix <- matrix(runif(d, -5, 5), ncol=d, nrow=1)
    theta <- theta_fix
    
    for (k in 0:K){
      if (k==1){theta <- matrix(theta[2,], nrow=1, ncol=d)}

      psi <- make_psi(theta)
       
        if (verbose && k %% 20 == 0) {
          e_xa <- e_n(A,X)
          mu_xa<-mu_n(A,X)
          nu_xa<- nu_n(A,X)
          psi_x <- psi(X)
          sb <- sigma_beta(psi,X, beta, centered)
          rho11<- (Y-mu_xa)*(2*A-1)
          rho22 <- (Xi -nu_xa)*(2*A-1)
            msg <- sprintf("FW-Est: iteration %i, value %f", k, J_funct_appprox(e_xa, psi_x, lambda, rho11, rho22, sb))
            message(msg)
        }
        theta_opt <- SGD_estimation(df, theta_fix, lambda, beta, centered, psi, e_n, mu_n, nu_n, (verbose && k %% 10 == 0))

        theta <- rbind(theta, theta_opt)
    }
    return(theta)
}


#' Debiasing procedure TMLE-like
#'
#' Debias the initial nuisance parameters estimation  
#' by calibrating \code{psi} such that the first term of the EIC of the objective function L 
#' becomes 0. Such calibration is performed via a 
#'
#' @param df A data frame n x d (input data).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to center the policy.
#' @param e_n The propensity score function, taking A and X as input.
#' @param mu_n The conditional expectation of Y function, taking A and X as input.
#' @param nu_n The conditional expectation of Xi function, taking A and X as input.
#' @param precision A numeric scalar that determines the convergence precision desired.
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the optimized parameter \code{theta},  
#'         where each row represents the k-th \code{theta} solution at iteration \code{k}.
#' @export
debias_procedure <- function(df, e.nj, mu.nj,nu.nj, lambda, beta,precision){
  X <- df %>% select(starts_with("X")) %>% as.matrix()

  res <- FW_estimation(df, lambda,  beta, centered, e.nj, mu.nj, nu.nj,precision)
  psi <- make_psi(res)

  Delta_mu_nj <- function(X){
    mu1_star <-mu.nj(1,X) + psi(X)
    mu0_star <- mu.nj(0,X) - psi(X)
    return(mu1_star-mu0_star)
  }

  Delta_nu_nj <- function(X){
    nu1_star <- nu.nj(1,X) + sigma_beta(psi,X, beta, centered)
    nu0_star <- nu.nj(0,X) - sigma_beta(psi,X, beta, centered)
    return(nu1_star-nu0_star)
  }
  return(list(Delta_mu_nj, Delta_nu_nj))
}