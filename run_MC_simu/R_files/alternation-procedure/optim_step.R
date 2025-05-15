args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])
fold <- as.integer(args[3])
k <- as.integer(args[4])  # iteration index

# Load necessary libraries and data
library(tidyverse)
folder <- readRDS("run_MC_simu/current_scenario.rds")
source(file.path("src", "algorithm_architecture.R"))
source(file.path(folder, "synthetic_data.R"))
source(file.path("src", "new_est_proced.R"))

df <- read.csv(file.path(folder, "results", "data", "estimated", paste0("df_", i, ".csv")))
X <- df %>% select(starts_with("X.")) %>% as.matrix()

beta <- readRDS(file.path(folder,"results","data","beta.rds"))
alpha <- 0.1
precision <- 0.025
lambda <- j
centered <- FALSE

mu.hat.nj <- readRDS(file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))
ps.hat <- readRDS(file.path(folder,"results","data","prop_score",paste0("ps.hat.nj_",i,".rds")))

mu0 <- mu.hat.nj[[fold]]
nu0 <- nu.hat.nj[[fold]]
prop_score <- ps.hat[[fold]]

# Define the previous and current step file paths
step_file_prev <- file.path(folder, "results", "new_estimated", paste0(i,j,fold,"_step_", k - 1, ".rds"))
step_file_current <- file.path(folder, "results", "new_estimated", paste0(i,j,fold,"_step_", k, ".rds"))

if (file.exists(step_file_prev)) {
  prev <- readRDS(step_file_prev)
} else {
  stop("Previous step file does not exist: ", step_file_prev)
}

# Extract variables from previous step
mu_XA <- prev$mu_XA
nu_XA <- prev$nu_XA
psi_vec <- as.matrix(prev$psi_vec)
sigma_psi_vec <- as.matrix(prev$sigma_psi_vec)
epsilon_1_vec <- prev$epsilon_1_vec
epsilon_2_vec <- prev$epsilon_2_vec

# Calculate the new epsilon values
H_cst <- HX(df$Treatment, X, prop_score)
new_epsilon_comb <- new_est(df$Y, df$Xi, mu_XA, nu_XA, psi_vec[, ncol(psi_vec)], sigma_psi_vec[, ncol(sigma_psi_vec)], H_cst)

epsilon_1_vec <- rbind(epsilon_1_vec, new_epsilon_comb[1])
epsilon_2_vec <- rbind(epsilon_2_vec, new_epsilon_comb[2])

# Update mu_XA and nu_XA
mu_XA <- df$Treatment * update_mu(1, X, mu0, epsilon_1_vec, psi_vec, prop_score) +
         (1 - df$Treatment) * update_mu(0, X, mu0, epsilon_1_vec, psi_vec, prop_score)
nu_XA <- df$Treatment * update_nu(1, X, nu0, epsilon_2_vec, sigma_psi_vec, prop_score) +
         (1 - df$Treatment) * update_nu(0, X, nu0, epsilon_2_vec, sigma_psi_vec, prop_score)

# Define new functions for delta_mu and delta_nu
Delta_mu <- function(X) { update_mu(1, X, mu0, epsilon_1_vec, psi_vec, prop_score) - update_mu(0, X, mu0, epsilon_1_vec, psi_vec, prop_score) }
Delta_nu <- function(X) { update_nu(1, X, nu0, epsilon_2_vec, sigma_psi_vec, prop_score) - update_nu(0, X, nu0, epsilon_2_vec, sigma_psi_vec, prop_score) }

# Calculate new theta and psi
theta <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision)
psi <- make_psi(theta)
psi_X <- psi(X)
psi_vec <- cbind(psi_vec, psi_X)

# Calculate the change between the previous and current psi_vec
psi_diff <- sum((psi_vec[,ncol(psi_vec)-1] - psi_X)^2)

# Calculate sigma and save results
sb <- sigma_beta(psi_X, beta, centered)
sigma_psi_vec <- cbind(sigma_psi_vec, sb)

out <- list(
  iter = k,
  mu_XA = mu_XA,
  nu_XA = nu_XA,
  psi_vec = psi_vec,
  sigma_psi_vec = sigma_psi_vec,
  epsilon_1_vec = epsilon_1_vec,
  epsilon_2_vec = epsilon_2_vec
)

# Save the new step
saveRDS(out, file = step_file_current)

# Delete the previous step file if it exists
if (file.exists(step_file_prev)) {
  file.remove(step_file_prev)
  cat("Deleted previous step file:", step_file_prev, "\n")
}

if (psi_diff < 0.25) {
  cat("Stopping iteration at k =", k, "because the difference is less than 0.5\n")
  cat("Outputting final theta:\n")
  print(theta)
  saveRDS(theta, file = file.path(folder, "results", "new_estimated", paste0(i,j,fold,"_final_theta.rds")))
  break
}