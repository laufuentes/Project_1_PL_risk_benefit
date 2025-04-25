args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
print(i)

folder <- readRDS("run_MC_simu/current_scenario.rds")
# Build filename
filename <- file.path(folder,"results","data","estimated",paste0("df_", i, ".csv"))


library(dplyr)
library(tidyverse)

df <- read.csv(filename, stringsAsFactors=FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()
n <- readRDS(file.path(folder,"results","data","n.rds"))
Jfold <-readRDS(file.path(folder,"results","data","Jfold.rds"))

source(file.path("src","estimation.R"))
source(file.path("src","cross-fitting.R"))

folds <- CVFolds(
    N = n,
    id = NULL,
    Y = df,
    cvControl = list(V = Jfold, stratifyCV = FALSE, shuffle = TRUE)
)

folds_df <- do.call(rbind, lapply(seq_along(folds), function(v){data.frame(index=folds[[v]], fold=v)}))
s <- folds_df$fold[order(folds_df$index)]
saveRDS(s, file = file.path(folder,"results","data","S",paste0("s_",i,".rds")))

initial_nparams <- nuissance_params(s,df)
mu.hat.nj <- initial_nparams[[2]]
saveRDS(mu.hat.nj, file = file.path(folder,"results","data","Mu",paste0("mu.hat.nj_",i,".rds")))

nu.hat.nj <- initial_nparams[[3]]
saveRDS(nu.hat.nj, file = file.path(folder,"results","data","Nu",paste0("nu.hat.nj_",i,".rds")))
