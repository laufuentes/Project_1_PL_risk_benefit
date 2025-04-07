library(tidyverse)
library(parallel)

set.seed("2025")
source("src/synthetic_data.R")
source("src/objective_functions.R")
source("src/algorithm_architecture.R")
source("src/optim_functions.R")
source("src/estimation.R")
source("src/plot_fcts.R")
expit <- plogis
logit <- qlogis

###########################
######## Parameters #######
###########################

n <- 2*1e3 # number of individuals 
alpha <- 0.1 # constraint tolerance
centered <- FALSE # centering of policy
epsilon <- 0.03 # early stopping parameter 
precision <- 0.025

#### Estimation parameters ####
Jfold<- 4 # number of K-folds
technique <- "probability.forest" # technique for estimation
cores_used<-8

#### Study simulation grid ####
B <- c(0.05,0.25,0.50,1,1.5,2) # beta candidates 
Lambda <- seq(0,10, 1)  # lambda candidates
param_combinations <- expand.grid(lambda = Lambda, beta = B)
param_combinations_est <- expand.grid(Fold=1:Jfold,lambda = Lambda, beta = B)

###########################
# Data generation process #
###########################
exp <- data_gen(n)
df_complete <- exp[[1]]
df <- exp[[2]]
X <- df %>% select(starts_with("X.")) %>% as.matrix()

synthetic_setting_plot(df_complete, delta_Y,delta_Xi)

##################################
## Oracular nuisance parameters ##
##################################
policies <- mclapply(1:nrow(param_combinations), function(i) {
    if(i%%100==0){print(i)}
  optimize_combination(i, param_combinations, delta_Y, delta_Xi)
}, mc.cores = cores_used-2, mc.preschedule=FALSE)

res_or<-mclapply(1:nrow(param_combinations), function(i) {
    process_policy(i,param_combinations,policies,X,
    list(df$y1,df$y0), delta_Y, delta_Xi,centered,alpha)
}, mc.cores = cores_used-2, mc.preschedule = FALSE) 

saveRDS(res_or, file = "opt_results/oracular/res.rds")

results_or <- as.data.frame(do.call(rbind, res_or))
write.csv(
  results_or %>% select(-optimal_x),
  paste0("opt_results/oracular/Results.csv")
)

idx_opt_obj <- which(
  results_or$obj == max(results_or$obj[results_or$constraint <= 0])
)
idx_opt <- idx_opt_obj

lambda_evol(
    results_or, 
    "oracular", 
    results_or$beta[[idx_opt]]
)

geom_points_fct(results_or,idx_opt, df, "oracular", centered)


#res <- readRDS("opt_results/oracular/res.rds")

##################################
# Nuissance parameter estimation #
##################################
s<- partition_data(Jfold, df)
## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
#e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[2]]
nu.hat.nj <- initial_nparams[[3]]

Delta_mu_nj_folds <- vector("list", Jfold)
Delta_nu_nj_folds <- vector("list", Jfold)
CC_mu <- vector("list", Jfold)
CC_nu<- vector("list", Jfold)

# Loop over each fold to compute the contrasts
for (fold in 1:Jfold) {
  # Get the mu and nu for the current fold
  mu.nj <- mu.hat.nj[[fold]]
  nu.nj <- nu.hat.nj[[fold]]
  
  # Define the Delta_mu and Delta_nu functions for the current fold
  Delta_mu_nj_folds[[fold]] <- function(X) { mu.nj(1, X) - mu.nj(0, X) }
  Delta_nu_nj_folds[[fold]] <- function(X) { nu.nj(1, X) - nu.nj(0, X) }
  CC_mu[[fold]]<- Delta_mu_nj_folds[[fold]](X)
  CC_nu[[fold]]<-Delta_nu_nj_folds[[fold]](X)
}

policies <- mclapply(1:nrow(param_combinations_est), function(i) {
    if(i%%100==0){print(i)}
  optimize_combination_Tlearner(i, param_combinations_est, Delta_mu_nj_folds, Delta_nu_nj_folds)
}, mc.cores = cores_used-2, mc.preschedule=FALSE)

res<-mclapply(1:nrow(param_combinations_est), function(i) {
    process_policy(i,param_combinations_est,policies,X,
    c(CC_mu[[param_combinations_est$Fold[[i]]]],
    CC_nu[[param_combinations_est$Fold[[i]]]]), 
    Delta_mu_nj_folds[[param_combinations_est$Fold[[i]]]], 
    Delta_nu_nj_folds[[param_combinations_est$Fold[[i]]]],
    centered,alpha)
}, mc.cores = cores_used-2, mc.preschedule = FALSE) 
saveRDS(res, file = "opt_results/estimation_T/res.rds")


results <- as.data.frame(do.call(rbind, res))
write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/estimation_T/Results.csv")
)

result_T <- results %>%
  group_by(lambda, beta) %>%
  summarise(
    obj = mean(obj),
    risk = mean(risk),
    constraint = mean(constraint),
    policy_value = mean(policy_value),
    .groups = "drop"
  )

write.csv(
  result_T,
  paste0("opt_results/estimation_T/Results.csv")
)

idx_opt_obj <- which(
  result_T$obj == max(result_T$obj[result_T$constraint <= 0])
)
idx_opt <- idx_opt_obj

lambda_evol(
    result_T, 
    "estimation_T", 
    result_T$beta[[idx_opt]]
)

theta_lambda0 <- Final_policy(0, 0.05,X, s, lambda<- 5
beta <- 0.55
theta_opt <- Final_policy(lambda, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds)
optimal_psi <- make_psi(theta_opt)
optimal_x <- optimal_psi(X))
psi0<- make_psi(theta_lambda0)
psi_x0 <- psi0(X)

lambda<- 5
beta <- 0.55
theta_opt <- Final_policy(lambda, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds)
optimal_psi <- make_psi(theta_opt)
optimal_x <- optimal_psi(X)

p0 <- gamma_plot_funct(psi_x0, 0, 0.05, df, centered)+ theme(legend.position = "none")
p <- gamma_plot_funct(optimal_x, lambda, 0.05, df, centered)+ theme(legend.position = "none")
opt_plots <- plot_grid(p0,p, ncol=2)
ggsave(paste0("images/estimation_T/2-Optimal_policy.pdf"),opt_plots, width = 8, height = 4)

sb_est <- sigma_beta(optimal_psi, X, 0.05, centered)
sb_or <- sigma_beta()

tibble(
  est=sb_est, 
  ora=sb_or
)%>%
ggplot()+
geom_point(aes(x=est,y=ora))+
geom_abline(intercept=0,slope=1, col="red")

ggsave("images/estimation/toto.pdf")