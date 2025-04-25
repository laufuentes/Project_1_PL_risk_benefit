source("src/synthetic_data.R")
source("src/cross-fitting.R")
library(dplyr)
library(ggplot2)
library(grf)
library(gridExtra)

n <- 1e4

df<- data_gen(n)[[2]]
X <- df %>% select(starts_with("X."))%>% as.matrix()
Jfold <-readRDS(file.path("MC","results","data","Jfold.rds"))

folds <- CVFolds(
    N = n,
    id = NULL,
    Y = df,
    cvControl = list(V = Jfold, stratifyCV = FALSE, shuffle = TRUE)
)

folds_df <- do.call(rbind, lapply(seq_along(folds), function(v){data.frame(index=folds[[v]], fold=v)}))
s <- folds_df$fold[order(folds_df$index)]


res_dir_mu_est <- file.path("MC", "results", "data", "Mu")
res_dir_nu_est <- file.path("MC", "results", "data", "Nu")

# List full file paths
res_files_mu_est <- list.files(path = res_dir_mu_est, full.names = TRUE)
res_files_nu_est <- list.files(path = res_dir_nu_est, full.names = TRUE)

# Extract and sort by lambda number for oracular files
numbers_mu <- as.numeric(sub(".*mu.hat.nj_(\\d+)\\.rds", "\\1", res_files_mu_est))
numbers_nu <- as.numeric(sub(".*nu.hat.nj_(\\d+)\\.rds", "\\1", res_files_nu_est))
res_files_mu_est <- res_files_mu_est[order(numbers_mu)]
res_files_nu_est <- res_files_nu_est[order(numbers_nu)]

all_mu <- lapply(res_files_mu_est, readRDS)
all_nu <- lapply(res_files_nu_est, readRDS)

CATE<- function(mu.hat.nj,X){
  Delta_mu_nj_folds <- lapply(mu.hat.nj, function(mu.nj) {
  function(X) mu.nj(1, X) - mu.nj(0, X)
})
  Delta_mu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      delta_mu <- Delta_mu_nj_folds[[fold]]
      out[s==fold] = delta_mu(X_fold)
    }
    return(out)
  }
  return(Delta_mu_CV)
}

CATE_Y_fct <- lapply(all_mu, function(mu)CATE(mu,X))
CATE_Y <- lapply(CATE_Y_fct, function(fct)fct(X))
CATE_df_Y <- do.call(cbind, lapply(CATE_Y, function(x) x))

CATE_Xi_fct <- lapply(all_nu, function(nu)CATE(nu,X))
CATE_Xi <- lapply(CATE_Xi_fct, function(fct)fct(X))
CATE_df_Xi <- do.call(cbind, lapply(CATE_Xi, function(x) x))

df_est <- cbind(df,CATE_Y=CATE_df_Y,CATE_Xi=CATE_df_Xi)%>%as.data.frame()
df_or <- cbind(df,CATE_Y=delta_mu(X),CATE_Xi =delta_nu(X))%>%as.data.frame()

df_std_CATE <- df_est %>%
  rowwise() %>%
  mutate(std_CATE_Y = sd(c_across(starts_with("CATE_Y"))))%>% 
  mutate(std_CATE_Xi = sd(c_across(starts_with("CATE_Xi"))))

p11<- ggplot(df_std_CATE, aes(x = X.1, y = X.2, color = std_CATE_Y)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c() +
  labs(title = "Delta Mu Variability Over Iterations",
       color = "Std Dev") +
  theme_minimal()

p12<- ggplot(df_std_CATE, aes(x = X.1, y = X.2, color = std_CATE_Xi)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c() +
  labs(title = "Delta Nu Variability Over Iterations",
       color = "Std Dev") +
  theme_minimal()

combined_plot <- grid.arrange(p11, p12, ncol =2)
ggsave(file.path("MC","figures","CATE_std.pdf"),combined_plot)


df_mean <- df_est %>%
  rowwise() %>%
  mutate(mean_CATE_Y = mean(c_across(starts_with("CATE_Y"))))%>% 
  mutate(mean_CATE_Xi = mean(c_across(starts_with("CATE_Xi"))))

p21<- ggplot(df_mean, aes(x = X.1, y = X.2, color = mean_CATE_Y)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(limits = c(-1, 1)) +
  labs(title = "Mean Delta Mu Over Iterations",
       color = "Mean Delta Mu") +
  theme_minimal()

p22<- ggplot(df_mean, aes(x = X.1, y = X.2, color = mean_CATE_Xi)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(limits = c(0, 1)) +
  labs(title = "Mean Delta Nu Over Iterations",
       color = "Mean Delta Nu") +
  theme_minimal()

p31<- ggplot(df_or, aes(x = X.1, y = X.2, color = CATE_Y)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(limits = c(-1, 1)) +
  labs(title = "Oracular Delta Mu",
       color = "Delta Mu") +
  theme_minimal()

p32<- ggplot(df_or, aes(x = X.1, y = X.2, color = CATE_Xi)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(limits = c(0, 1)) +
  labs(title = "Oracular Delta Nu",
       color = "Delta Nu") +
  theme_minimal()

combined_plot_mean <- grid.arrange(p21, p22, p31,p32, nrow =2)
ggsave(file.path("MC","figures","CATE_mean.pdf"),combined_plot_mean)

########################
###### Densities #######
########################
oracular_Y_vec <- as.vector(delta_mu(X))
estimated_Y_vec <- as.vector(CATE_df_Y)

densities_Y <- data.frame(
  value = c(oracular_Y_vec, estimated_Y_vec),
  source = c(
    rep("Oracular", length(oracular_Y_vec)),
    rep("Estimated", length(estimated_Y_vec))
  )
)

pY <- ggplot(densities_Y, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density of Delta Mu comparison",
       x = "Function Output",
       y = "Density") +
  theme_minimal()

oracular_Xi_vec <- as.vector(delta_nu(X))
estimated_Xi_vec <- as.vector(CATE_df_Xi)

densities_Xi <- data.frame(
  value = c(oracular_Xi_vec, estimated_Xi_vec),
  source = c(
    rep("Oracular", length(oracular_Xi_vec)),
    rep("Estimated", length(estimated_Xi_vec))
  )
)

# Plot the density
pXi <- ggplot(densities_Xi, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density of Delta Nu comparison",
       x = "Function Output",
       y = "Density") +
  theme_minimal()

combined_plot_density <- grid.arrange(pY, pXi, ncol =2)
ggsave(file.path("MC","figures","density_comparison.pdf"),combined_plot_density)

