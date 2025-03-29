library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(madness)

source("src/synthetic_data.R")
source("src/tool_box.R")
source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/plot_fcts_change_after.R")

expit <- plogis
logit <- qlogis

##############
# Parameters #
##############

# GENERAL PARAMETERS
n <- 1e3 # number of individuals 
alpha <- 0.1 # constraint tolerance
centered <- TRUE # centering of policy
epsilon <- 0.03 # early stopping parameter 
precision <- 0.025

# DATA GENERATION PARAMETERS
setting <- "Other_1" # setting selected
option <- option_det(setting, "_") 

# GRID SEARCH PARAMETERS
B <- seq(0.05, 2, 0.05) # beta candidates 
Lambda <- seq(0,10, 1)  # lambda candidates

param_combinations <- expand.grid(lambda = Lambda, beta = B)

###########################
# Data generating process #
###########################
exp <- data_gen(n, option)
df <- exp[[1]]

df$Xi_0 <- ifelse(df$p0>2,1,0)
df$Xi_1 <- ifelse(df$p1>2,1,0)
covariates <- df %>% select(starts_with("X")) %>% as.matrix()

##########################
# Optimization functions #
##########################
optimize_combination <- function(i, param_combinations){
  thetas <- FW(covariates, 
  param_combinations$lambda[i], 
  param_combinations$beta[i], 
  alpha, delta_Y, delta_Xi, precision)

  return(thetas) 
}

policies <- mclapply(1:nrow(param_combinations), function(i) {
    if(i%%100==0){print(i)}
  optimize_combination(i, param_combinations)
}, mc.cores = detectCores(),mc.preschedule = FALSE)

process_policy <- function(
    idx,
    param_combinations,
    thetas,
    covariates,
    df_complete,
    centered,
    alpha) {
    # Extract the policy for the current index
    theta <- thetas[[idx]]
    psi <-make_psi(theta)
    results <- data.frame(
        lambda = param_combinations$lambda[idx],
        beta = param_combinations$beta[idx],
        optimal_x = I(list(psi(covariates))), # I() wraps the list to avoid issues with data frames
        risk = R_p(psi, covariates, delta_Y),
        constraint = S_p(
            psi,
            covariates,
            param_combinations$beta[idx],
            alpha, centered, delta_Xi
            ),
        obj = L(psi, covariates,param_combinations$lambda[idx], param_combinations$beta[idx], alpha, centered, delta_Y, delta_Xi),
        policy_value = policy_values(
          psi,
          covariates, 
          c(df_complete$y1, df_complete$y0),
          param_combinations$beta[idx],
          centered,
          alpha
        )
    )
    colnames(results) <- c(
        "lambda",
        "beta",
        "optimal_x",
        "risk",
        "constraint",
        "obj",
        "policy_value")
    return(results) # Return the updated results for this index
}


save_policies_to_csv <- function(policies, file_name = "opt_results/oracular/policies_output.csv") {
  if (is.list(policies)) {
    policies_df <- do.call(rbind, lapply(policies, as.data.frame))  # Convert list to data frame
    write.csv(policies_df, file = file_name, row.names = FALSE)
    message("Policies saved to ", file_name)
  } else {
    stop("Policies output is not a list and cannot be converted to a data frame.")
  }
}

save_policies_to_csv(policies)
saveRDS(policies, file = "opt_results/oracular/Thetas.rds")


res<-mclapply(1:nrow(param_combinations), function(i) {
    process_policy(i,param_combinations,policies,covariates,
    df,centered,alpha)
}, mc.cores = detectCores(), mc.preschedule = FALSE) 

saveRDS(res, file = "opt_results/oracular/res.rds")
results <- as.data.frame(do.call(rbind, res))

# Save the results to a CSV file
write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/oracular/", setting, ".csv")
)


idx_opt_pol <- which(
  results$policy_value == max(results$policy_value[results$constraint <= 0])
)
idx_opt_obj <- which(results$obj == max(results$obj[results$constraint <= 0]))


idx_opt <- idx_opt_pol[1]

lambda_evol(
  results %>% select(-optimal_x),
  "oracular",option, 
  results$lambda[idx_opt],
  results$beta[idx_opt])

geom_points_fct(results, idx_opt, df, "oracular", option, centered=TRUE)

gamma_lambda_plot(results, df, "oracular", option)

##################################
######### Plot functions #########
##################################



sigma_beta_fixed_X <- function(psi_X){
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered) {
        cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi_X)) / (1 + exp(-beta))) + cent
}

L_fixed_X <-function(psi_X){
  mean(psi_X^2 -2*psi_X*delta_Y(X) + lambda*(sigma_beta_fixed_X(psi_X)*delta_Z(X)-alpha))
}


theta <- matrix(runif(ncol(X),-1,1),nrow=1, ncol=ncol(X))
psi_opt_X <- optim(L_fixed_X, par=rep(0,n), method="L-BFGS-B", lower=-1, upper=1)$par

policy_free <- sigma_beta_fixed_X(psi_opt_X)

tibble(
  FW=policy_FW,
  free=policy_free
) %>%
  ggplot() +
  geom_point(aes(x=FW, y=free))+geom_abline(intercept=0,slope=1, color="red", linetype="dashed")

