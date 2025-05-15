library(ggplot2)
library(tidyverse)
library(grf)

folder <- readRDS("run_MC_simu/current_scenario.rds")

source("src/optim_functions.R")
source(file.path(folder,"synthetic_data.R"))

i<- 1
l<-1
fold<-2

df <- read.csv(file.path(folder,"results","data","estimated",paste0("df_",i,".csv")))
X <- df %>% select(starts_with("X."))%>% as.matrix()

theta_opt <- readRDS(file.path(folder,"results","oracular","theta_opt","opt",paste0("theta_opt_",i,".rds")))
psi_opt <- make_psi(theta_opt)
psi_opt_X <- psi_opt(X)
sigma_psi_opt <- sigma_beta(psi_opt_X)

theta_new <- readRDS(file.path(folder,"results","new_estimated",paste0(i,l,fold,"_final_theta.rds")))
psi_new <- make_psi(theta_new)
psi_new_X <- psi_new(X)
sigma_psi_new <- sigma_beta(psi_new_X)


theta_tlr <- readRDS(file.path(folder,"results","estimated","theta_opt",paste0("theta_opt_",i,".rds")))
psi_tlr <- make_psi(theta_tlr)
psi_tlr_X <- psi_tlr(X)
sigma_psi_tlr <- sigma_beta(psi_tlr_X)


df <- data.frame(opt=psi_opt_X, tlr=psi_tlr_X, corrected=psi_new_X)
df_long <- bind_rows(
  data.frame(x = df$opt, y = df$tlr, pair = "opt vs tlr"),
  data.frame(x = df$opt, y = df$corrected, pair = "opt vs tlr-corrected")
)
comp<- ggplot(df_long, aes(x = x, y = y, color = pair)) +
  geom_point(alpha = 0.7) + geom_abline(intercept=0,slope=1,color="red")+
  labs(x = "opt", y = "tlr or tlr-corrected", title = "Comparison of opt with tlr and tlr-corrected") +
  theme_minimal()

ggsave(comp,filename=file.path(folder,"figures","new_estimated",paste0("Comparison_",i,l,fold,".pdf")))


prefix<- paste0(i,l,fold, "_step_")
target_dir <- file.path(folder,"results","new_estimated")
file_path<- list.files(pattern = paste0("^", prefix, "[0-9]+\\.rds$"), path = target_dir)
res <- readRDS(file.path(target_dir,file_path))
diff_sums <- colSums((res$psi_vec[,ncol(res$psi_vec)] - res$psi_vec[,-1])^2)

diff_plot <- ggplot()+geom_point(aes(x=1:length(diff_sums),y=diff_sums))
ggsave(diff_plot,filename=file.path(folder,"figures","new_estimated",paste0("Loss_",i,l,fold,".pdf")))