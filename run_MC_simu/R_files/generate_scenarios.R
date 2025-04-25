args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]

saveRDS(folder, file = "run_MC_simu/current_scenario.rds")
expit <- plogis
logit <- qlogis
library(parallel)
library(tidyverse)

# Load sources 
source(file.path(folder,"synthetic_data.R"))


n <- 2*1e3
saveRDS(n, file = file.path(folder,"results","data","n.rds"))

MC_iterations <- 1e2
saveRDS(MC_iterations, file = file.path(folder,"results","data","MC_iter.rds"))

for (b in 1:MC_iterations){
    exp <- data_gen(n,seed=1e3+b)
    df_complete <- exp[[1]]
    write.csv(df_complete,file.path(folder,"results","data","oracular",paste0("df_complete_", b,".csv")))
    df <- exp[[2]]
    write.csv(df,file.path(folder,"results","data","estimated",paste0("df_",b,".csv")))
}

source("src/plot_functions.R")
synthetic_setting_plot(df_complete, delta_mu, delta_nu, paste0(file.path(folder,"figures")))

Lambda <- seq(0,12,1)
saveRDS(Lambda, file = file.path(folder,"results","data","Lambda.rds"))

beta <- 0.05
saveRDS(beta, file = file.path(folder,"results","data","beta.rds"))

Jfold <- 5
saveRDS(Jfold, file = file.path(folder,"results","data","Jfold.rds"))

param_combinations <- expand.grid(lambda=Lambda, beta = beta)
saveRDS(param_combinations, file = file.path(folder,"results","data","param_combinations.rds"))