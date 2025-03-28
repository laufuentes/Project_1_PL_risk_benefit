setwd("~/Project_1_PL_risk_benefit/")
library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(numDeriv)

source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/synthetic_data.R")
source("src/estimators.R")
source("src/tool_box.R")
source("src/plot_fcts_change_after.R")
source("src/algorithm_architecture.R")

######################
# General parameters #
######################
n <- 1e3 # number of individuals 
alpha <- 0.1 # constraint tolerance
centered <- TRUE # centering of policy
epsilon <- 0.03 # early stopping parameter 
precision <- 0.025

##  Data generation parameters
setting <- "Other_1" # setting selected
option <- option_det(setting, "_") 

# Grid search parameters
#B <- c(0.05,0.1)
#Lambda <- c(0,5,10)
B <- seq(0.05, 2, 0.05) # beta candidates 
Lambda <- seq(0,15, 1) # lambda candidates

# Estimation parameters
Jfold<- 5 # number of K-folds
technique <- "cv.glmnet" # technique for estimation

## Data generating process ##
exp <- data_gen(n, option)
df <- exp[[2]]
df$Xi <- ifelse(df$Z>2,1,0)

X <- df %>% select(starts_with("X")) %>% as.matrix()

param_combinations <- expand.grid(Fold=1:Jfold,lambda = Lambda, beta = B)

# Step 1: Obtain initial estimates of nuisance parameters
## Partition dataset 
s<- partition_data(Jfold, df)
## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[3]]
nu.hat.nj <- initial_nparams[[4]]

combination_results <- mclapply(1:nrow(param_combinations), function(i) {
    if(i%%100==0){print(i)}
  algo_combination(param_combinations[i,])
}, mc.cores = detectCores(),mc.preschedule = FALSE)


RowCountsDF <- data.frame(
  CombinationID = 1:length(combination_results),  # Unique ID for each combination
  NumRows = sapply(combination_results, function(res) nrow(res[[2]])) # Extract row count
)

Theta_df <- as.data.frame(do.call(rbind, lapply(combination_results, `[[`, 2)))

save_policies_to_csv(Theta_df)

Thetas <- lapply(combination_results, `[[`, 2)
Causal_contrast <- lapply(combination_results, `[[`, 1)

# Save Thetas
saveRDS(Thetas, file = "opt_results/Thetas.rds")
# Save Causal_contrast
saveRDS(Causal_contrast, file = "opt_results/Causal_contrast.rds")

res<-mclapply(1:nrow(param_combinations), function(i) {
    process_policy(i,param_combinations,Thetas,X,mu.hat.nj,
    Causal_contrast,
    centered,alpha)
}, mc.cores = detectCores(), mc.preschedule = FALSE) 

results <- as.data.frame(do.call(rbind, res))

# Save the results to a CSV file
write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/estimation/", setting, ".csv")
)


idx_opt_pol <- which(
  results$policy_value == max(results$policy_value[results$constraint <= 0])
)
idx_opt_obj <- which(results$obj == max(results$obj[results$constraint <= 0]))


idx_opt <- idx_opt_pol[1]

lambda_evol(
  results %>% filter(fold==results$fold[idx_opt])%>% select(-optimal_x),
  "estimation",option, 
  results$lambda[idx_opt],
  results$beta[idx_opt])


optimal_x <- as.numeric(results[idx_opt, ]$optimal_x[[1]])
p <- gamma_plot_funct(optimal_x, results$lambda[idx_opt], results$beta[idx_opt], df, option, centered)
ggsave("images/estimation/optimal_policy.pdf",p)

#geom_points_fct(results%>%filter(fold==results$fold[[idx_opt]]), idx_opt, df, "estimation", option, centered=TRUE)

#gamma_lambda_plot(results, df, "estimation", option)


# To load them back later:
# loaded_Thetas <- readRDS("Thetas.rds")
# loaded_Causal_contrast <- readRDS("Causal_contrast.rds")