args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
print(i)
# Build filename
filename <- file.path("MC","results","data","estimated",paste0("df_", i, ".csv"))

library(dplyr)
library(tidyverse)

df <- read.csv(filename, stringsAsFactors=FALSE)
X <- df %>% select(starts_with("X.")) %>% as.matrix()
n <- readRDS(file.path("MC","results","data","n.rds"))
Jfold <-readRDS(file.path("MC","results","data","Jfold.rds"))

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
saveRDS(s, file = file.path("MC","results","data","S",paste0("s_",i,".rds")))

initial_nparams <- nuissance_params(s,df)
mu.hat.nj <- initial_nparams[[2]]
saveRDS(mu.hat.nj, file = file.path("MC","results","data","Mu",paste0("mu.hat.nj_",i,".rds")))

nu.hat.nj <- initial_nparams[[3]]
saveRDS(nu.hat.nj, file = file.path("MC","results","data","Nu",paste0("nu.hat.nj_",i,".rds")))


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

saveRDS(Delta_mu_nj_folds, file = file.path("MC","results","data","Delta_Mu",paste0("Delta_mu_nj_folds_",i,".rds")))
saveRDS(Delta_nu_nj_folds, file = file.path("MC","results","data","Delta_Nu",paste0("Delta_nu_nj_folds_",i,".rds")))
saveRDS(CC_mu, file = file.path("MC","results","data","CC_mu",paste0("CC_mu_",i,".rds")))
saveRDS(CC_nu, file = file.path("MC","results","data","CC_nu",paste0("CC_nu_",i,".rds")))