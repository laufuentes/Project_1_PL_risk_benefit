args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])
fold <- as.integer(args[3])

# Load packages and data
library(tidyverse)
folder <- readRDS("run_MC_simu/current_scenario.rds")
source(file.path("src", "algorithm_architecture.R"))
source(file.path(folder, "synthetic_data.R"))
source(file.path("src", "new_est_proced.R"))

df <- read.csv(file.path(folder, "results", "data", "estimated", paste0("df_", i, ".csv")))
X <- df %>% select(starts_with("X.")) %>% as.matrix()
s <- readRDS(file.path(folder,"results","data","S",paste0("s_",i,".rds")))
mu.hat.nj <- readRDS(file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))
nu.hat.nj <- readRDS(file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))
ps.hat <- readRDS(file.path(folder,"results","data","prop_score",paste0("ps.hat.nj_",i,".rds")))
beta <- readRDS(file.path(folder,"results","data","beta.rds"))
alpha <- 0.1
precision <- 0.05
lambda <- j
centered <- FALSE

mu0 <- mu.hat.nj[[fold]]
nu0 <- nu.hat.nj[[fold]]
prop_score <- ps.hat[[fold]]

Delta_mu <- function(X) { mu0(1, X) - mu0(0, X) }
Delta_nu <- function(X) { nu0(1, X) - nu0(0, X) }

theta_init <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision)
psi <- make_psi(theta_init)
psi_X <- psi(X)
sb <- sigma_beta(psi_X, beta, centered)

out_list <- list(
  iter = 0,
  mu_XA = df$Treatment * mu0(1, X) + (1 - df$Treatment) * mu0(0, X),
  nu_XA = df$Treatment * nu0(1, X) + (1 - df$Treatment) * nu0(0, X),
  psi_vec = psi_X,
  sigma_psi_vec = sb,
  epsilon_1_vec = NULL,
  epsilon_2_vec = NULL
)

saveRDS(out_list, file = file.path(folder, "results", "new_estimated", paste0(i,j,fold,"_step_0.rds")))