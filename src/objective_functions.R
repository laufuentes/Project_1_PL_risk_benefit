
#' Probability of Treating Function
#'
#' Computes the probability of treatment assignment based on the input function \code{psi(X)},  
#' a parameter \code{beta}, and an optional centering adjustment.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.  
#'        Typically, \code{psi(X)} is a convex combination of past solutions in the Frank-Wolfe algorithm.
#' @param X A numeric matrix of size n x d (input data).
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering to the probability computation.
#'
#' @return A numeric vector of probabilities, where each element represents the treatment probability  
#'         for a corresponding row in \code{X}.
#' @export
sigma_beta <- function(psi, X, beta, centered) {
    c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered) {
        cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi(X))) / (1 + exp(-beta))) + cent
return(out)
}

#' Derivative of Probability of Treating Function
#'
#' Computes the derivative of the probability of treatment function \code{sigma_beta},  
#' which represents the sensitivity of treatment probability to changes in \code{psi} at X.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#'
#' @return A numeric vector representing the derivative of \code{sigma_beta} with respect to \code{psi} at X.
#' @export
sigma_beta_prime <- function(psi, X, beta){
    if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
    c_beta <- 1 / log(
      (1 + exp(beta)) / (1 + exp(-beta))
      )
    out <- c_beta *(beta*exp(beta*psi(X)))/(1+ exp(beta*psi(X)))
    return(out)
}

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
#' @param delta_Xi A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#'
#' @return A numeric scalar representing the constraint function value.
#' @export
S_p <- function(psi, X, beta, alpha, centered, delta_Xi){
    out <- mean(sigma_beta(psi,X, beta, centered) * delta_Xi(X)) - alpha
    return(out)
}

#' Objective Function
#'
#' Computes the objective function \eqn{L}, which balances the risk function \eqn{R_p}  
#' and the constraint function \eqn{S_p} using a parameter \code{lambda}.
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param alpha A numeric scalar representing the constraint tolerance.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#' @param delta_Y A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Xi A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#'
#' @return A numeric scalar representing the objective function value.
#' @export
L <- function(psi, X,lambda, beta, alpha, centered, delta_Y, delta_Xi){
    out <- R_p(psi, X,delta_Y) + lambda*S_p(psi, X, beta, alpha, centered, delta_Xi)
    return(out)
}

#' Gradient of the Objective Function
#'
#' Computes the gradient of the objective function \code{L} with respect to \code{psi} at X.  
#' The gradient is used in optimization algorithms like Stochastic Gradient Descent (SGD).
#'
#' @param psi A function that takes a numeric matrix \code{X} as input and returns a numeric vector.
#' @param X A numeric matrix of size n x d (input data).
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#' @param delta_Y A function of \code{X} that determines the difference between primary counterfactual outcomes.
#' @param delta_Xi A function of \code{X} that determines the difference between secondary counterfactual outcomes.
#'
#' @return A numeric vector representing the gradient of the objective function with respect to \code{psi(X)}.
#' @export
gradL <- function(psi, X, lambda, beta, centered, delta_Y, delta_Xi){
  2*(psi(X) - delta_Y(X)) + lambda*sigma_beta_prime(psi,X, beta)*delta_Xi(X)
}

#' Expected Policy Outcome
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
#'
#' @return A numeric scalar representing the expected outcome under the policy.
#' @export
policy_values <- function(psi, X,counterfactual_outcomes,beta,centered,alpha){
  sigma_psi <-sigma_beta(psi, X, beta, centered)
  policy <- rbinom(nrow(X), 1, sigma_psi)
  y1 <- counterfactual_outcomes[1]
  y0 <- counterfactual_outcomes[2]
  out <- mean(policy * y1 + (1 - policy) * y0)
  return(out)
}