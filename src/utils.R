#' Link Function
#'
#' Link function mapping \eqn{[-1,1]} to \eqn{[0,1]}, parametrized    
#' by \code{beta} with an optional centering.
#'
#' @param t A vector of numerics.
#' @param beta A numeric scalar (0.05 by default) controlling the sharpness of the link function.
#' @param centered A logical (FALSE by default) indicating whether to apply centering to the link function.
#'
#' @return A numeric vector of probabilities.
#' @export
sigma_beta <- function(t, beta=0.05, centered=FALSE) {
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
  if (centered) {
    cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
  } else {
    cent <- 0
  }
  out <- c_beta * log((1 + exp(beta * t)) / (1 + exp(-beta))) + cent
  return(out)
}

#' Derivative of Link Function
#'
#' Computes the derivative of the link function \code{sigma_beta},  
#' with respect to t.
#'
#' @param t A vector of numerics.
#' @param beta A numeric scalar (0.05 by default) controlling the sharpness of the link function.
#' @param centered A logical (FALSE by default) indicating whether to apply centering to the link function.
#'
#' @return The derivative of \code{sigma_beta} evaluated at t.
#' @export
sigma_beta_prime <- function(t, beta=0.05, centered=FALSE){
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  c_beta <- 1 / log(
    (1 + exp(beta)) / (1 + exp(-beta))
  )
  out <- c_beta *(beta*exp(beta*t))/(1+ exp(beta*t))
  return(out)
}


#' Oracular Approximation of Value Function
#'
#' Computes the expected outcome under a policy determined by the previously optimized \code{psi(X)}.  
#' The policy assigns treatment probabilistically based on \code{sigma_beta(psi(X))},  
#' and the expected outcome is calculated using counterfactual outcomes.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param counterfactual_outcomes A list of length 2 containing:
#'   \itemize{
#'     \item \code{counterfactual_outcomes[[1]]}: Outcome if treated (\eqn{Y(1)}).
#'     \item \code{counterfactual_outcomes[[2]]}: Outcome if not treated (\eqn{Y(0)}).
#'   }
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#' @param alpha A numeric scalar representing the constraint tolerance (not directly used here but for consistency).
#' @param B Integer (default 1e4) number of Monte Carlo repetitions.
#' @param seed Integer or NA (default value).
#'
#' @return A numeric scalar representing the expected outcome under the policy.
#' @export
V_p <- function(psi, beta, centered, alpha, B=1e4, seed=NA){
  `%>%`<- magrittr::`%>%`
  df <- data_gen(B,seed=seed)[[1]]
  X <- df%>%dplyr::select(dplyr::starts_with("X."))%>% as.matrix()
  y1 <- df$y1
  y0 <- df$y0
  
  psi_X <- psi(X)
  sigma_psi <-sigma_beta(psi_X, beta, centered)
  if(!is.na(seed)){
    set.seed(seed)
  }
  action <- stats::rbinom(nrow(X), 1, sigma_psi)
  out <- mean(action * y1 + (1 - action) * y0)
  return(out)
}

V_n <- function(psi,beta,mu.nj,alpha,B=1e4, seed=NA){
  `%>%`<- magrittr::`%>%`
  df <- data_gen(B,seed=seed)[[2]]
  X <- df%>%dplyr::select(dplyr::starts_with("X."))%>% as.matrix()
  y1 <- mu.nj(1,X)
  y0 <- mu.nj(0,X)
  
  psi_X <- psi(X)
  sigma_psi <-sigma_beta(psi_X, beta, centered)
  if(!is.na(seed)){
    set.seed(seed)
  }
  action <- stats::rbinom(nrow(X), 1, sigma_psi)
  out <- mean(action * y1 + (1 - action) * y0)
  return(out)
}