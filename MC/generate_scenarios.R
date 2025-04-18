library(parallel)
library(tidyverse)

# Load sources 
source(file.path("src","synthetic_data.R"))

n <- 2*1e3
saveRDS(n, file = file.path("MC","results","data","n.rds"))

MC_iterations <- 1e2
saveRDS(MC_iterations, file = file.path("MC","results","data","MC_iter.rds"))

for (b in 1:MC_iterations){
    exp <- data_gen(n,seed=1e3+b)
    df_complete <- exp[[1]]
    write.csv(df_complete,file.path("MC","results","data","oracular",paste0("df_complete_", b,".csv")))
    df <- exp[[2]]
    write.csv(df,file.path("MC","results","data","estimated",paste0("df_",b,".csv")))
}

Lambda <- seq(0,8,0.5)
saveRDS(Lambda, file = file.path("MC","results","data","Lambda.rds"))

beta <- 0.05
saveRDS(beta, file = file.path("MC","results","data","beta.rds"))

Jfold <- 5
saveRDS(Jfold, file = file.path("MC","results","data","Jfold.rds"))

param_combinations <- expand.grid(lambda=Lambda, beta = beta)
saveRDS(param_combinations, file = file.path("MC","results","data","param_combinations.rds"))
