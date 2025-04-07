source("src/synthetic_data.R")
expit <- plogis
logit <- qlogis

n <- 2*1e3
exp <- data_gen(n)
df_complete <- exp[[1]]
df <- exp[[2]]
X <- df %>% select(starts_with("X.")) %>% as.matrix()


write.csv(df_complete,paste0("opt_results/oracular/df_complete.csv"))
write.csv(df,paste0("opt_results/estimation_T/df.csv"))


lambda <- seq(0,10,1) 
saveRDS(lambda, file = "opt_results/data/Lambda.rds")

B <- c(0.05, 0.25,0.5,1,1.5,2)
saveRDS(B, file = "opt_results/data/B.rds")

Jfold<- 5
saveRDS(Jfold, file = "opt_results/data/Jfolds.rds")


source("src/estimation.R")
s<- partition_data(Jfold, df)
saveRDS(s, file = "opt_results/data/s.rds")

## Construct cross-fitting estimates
initial_nparams <- nuissance_params(s,df, technique)
#e.hat.nj <- initial_nparams[[2]]
mu.hat.nj <- initial_nparams[[2]]
saveRDS(mu.hat.nj, file = "opt_results/data/mu.hat.nj.rds")

nu.hat.nj <- initial_nparams[[3]]
saveRDS(nu.hat.nj, file = "opt_results/data/nu.hat.nj.rds")

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

saveRDS(Delta_mu_nj_folds, file = "opt_results/data/Delta_mu_nj_folds.rds")
saveRDS(Delta_nu_nj_folds, file = "opt_results/data/Delta_nu_nj_folds.rds")
saveRDS(CC_mu, file = "opt_results/data/CC_mu.rds")
saveRDS(CC_nu, file = "opt_results/data/CC_nu.rds")
