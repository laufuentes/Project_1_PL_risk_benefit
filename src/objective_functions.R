#' Risk Function for Conditional Average Treatment Effect (CATE)
#'
#' Computes the risk function \eqn{R_p} for estimating the Conditional Average Treatment Effect (CATE).  
#' The function minimizes the squared error between \code{psi(X)} and \code{delta_Y(X)}.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param delta_Y A function of \code{X} that determines the difference between primary counterfactual outcomes.
#'
#' @return A numeric scalar representing the risk function value.
#' @export
R_p <- function(psi, X, delta_Y){
    out <- mean(psi(X)^2 - 2 * psi(X)* delta_Y(X))
    return(out)
}

#' Constraint Function
#'
#' Computes the constraint function \eqn{S_p}, which ensures that the learned policy satisfies a constraint. 
#' This function enforces a limit on the expected impact of treatment via \code{delta_Z}.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param alpha A numeric scalar representing the constraint tolerance.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#' @param delta_Nu A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#'
#' @return A numeric scalar representing the constraint function value.
#' @export
S_p <- function(psi, X, beta, alpha, centered, delta_Nu){
  psi_X <- psi(X)
  out <- mean(sigma_beta(psi_X, beta, centered) * delta_Nu(X)) - alpha
    return(out)
}

#' Objective Function taking the form of a Lagrangian
#'
#' Computes the objective function, which balances the risk function \code{R_p}  
#' and the constraint function \code{S_p} using a parameter \code{lambda}.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Nu A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param alpha A numeric scalar (0.05 by default) representing the constraint tolerance.
#' @param beta A numeric scalar (0.05 by default) controlling the sharpness of the probability function.
#' @param centered A logical (FALSE by default) indicating whether to apply centering in \code{sigma_beta}.
#'
#' @return A numeric scalar representing the objective function value.
#' @export
Lagrangian_p <- function(psi, X, delta_Mu, delta_Nu, lambda, alpha=0.05, beta=0.05, centered=FALSE){
    out <- R_p(psi, X,delta_Mu) + lambda*S_p(psi, X, beta, alpha, centered, delta_Nu)
    return(out)
}

#' Gradient of the Objective Function
#'
#' Computes the gradient of the objective function with respect to \code{psi} at X.  
#' The gradient is used in optimization algorithms like Stochastic Gradient Descent (SGD).
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Nu A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param alpha A numeric scalar (0.05 by default) representing the constraint tolerance.
#' @param beta A numeric scalar (0.05 by default) controlling the sharpness of the probability function.
#' @param centered A logical (FALSE by default) indicating whether to apply centering in \code{sigma_beta}.
#'
#' @return A numeric vector representing the gradient of the objective function with respect to \code{psi(X)}.
#' @export
grad_Lagrangian_p <- function(psi, X, delta_Mu, delta_Nu, lambda, alpha=0.05, beta=0.05, centered=FALSE){
  psi_X <- psi(X)
  2*(psi(X) - delta_Mu(X)) + lambda*sigma_beta_prime(psi_X, beta, centered)*delta_Nu(X)
}