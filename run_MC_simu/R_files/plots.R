library(tidyverse)
library(ggplot2)
library(tidyr)

folder <- readRDS("run_MC_simu/current_scenario.rds")
########################
# Lambda evol #
########################

# Set up directories
res_dir_or <- file.path(folder, "results", "oracular", "group_by_lambda")
res_dir_est <- file.path(folder, "results", "estimated", "oracular", "group_by_lambda")

# List full file paths
res_files_or <- list.files(path = res_dir_or, full.names = TRUE)
res_files_est <- list.files(path = res_dir_est, full.names = TRUE)

# Extract and sort by lambda number for oracular files
lambda_nums_or <- as.numeric(sub(".*lambda_(\\d+)\\.csv", "\\1", res_files_or))
res_files_or <- res_files_or[order(lambda_nums_or)]
all_results_or <- lapply(res_files_or, function(x)read.csv(x, stringsAsFactors=FALSE))
results_or <- do.call(rbind, all_results_or)

# Extract and sort by lambda number for estimated files
lambda_nums_est <- as.numeric(sub(".*lambda_(\\d+)\\.csv", "\\1", res_files_est))
res_files_est <- res_files_est[order(lambda_nums_est)]
all_results_est <- lapply(res_files_est, function(x)read.csv(x, stringsAsFactors=FALSE))
results_est<- do.call(rbind, all_results_est)

# Reshape the data to long format (necessary for ggplot)
results_or_long <- results_or %>%
  gather(key = "variable", value = "value", policy_value, obj, risk, constraint)

results_est_long <- results_est %>%
  gather(key = "variable", value = "value", policy_value, obj, risk, constraint)

# Create separate boxplots for each variable, with color fill based on the variable
l_evol_or<- ggplot(results_or_long, aes(x = factor(lambda), y = value, fill = variable)) + 
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +  # Separate boxplots for each variable
  labs(x = "Lambda", y = "Value", title = "Boxplots of Policy Value, Objective, Risk, and Constraint by Lambda") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Choose a color palette
ggsave(paste0(folder,"/figures/oracular/lambda_evol_or.pdf"), l_evol_or)

l_evol_est<- ggplot(results_est_long, aes(x = factor(lambda), y = value, fill = variable)) + 
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +  # Separate boxplots for each variable
  labs(x = "Lambda", y = "Value", title = "Boxplots of Policy Value, Objective, Risk, and Constraint by Lambda") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Choose a color palette
ggsave(paste0(folder,"/figures/estimated/lambda_evol_est.pdf"), l_evol_est)


# Assuming results_or and results_est are already loaded into your environment

# Add a new column to each dataset to distinguish between 'or' and 'est'
results_or$group <- 'or'
results_est$group <- 'est'

# Combine both datasets into one
results_combined <- bind_rows(results_or, results_est)

# Reshape the data to long format (necessary for ggplot)
results_combined_long <- results_combined %>%
  gather(key = "variable", value = "value", policy_value, obj, risk, constraint)

# Create combined boxplot with color distinction between 'or' and 'est'
l_evol_both <- ggplot(results_combined_long, aes(x = factor(lambda), y = value, fill = group)) + 
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +  # Separate boxplots for each variable
  labs(x = "Lambda", y = "Value", title = "Boxplots of Policy Value, Objective, Risk, and Constraint by Lambda") +
  theme_minimal() +
  scale_fill_manual(values = c("or" = "blue", "est" = "red"))  # Manually set colors for 'or' and 'est'

ggsave(paste0(folder,"/figures/lambda_evol_both.pdf"), l_evol_both)

########################
# Psi opt #
########################
source("src/optim_functions.R")
source("src/utils.R")
source(file.path(folder,"synthetic_data.R"))

n <- 1e4
df <- data_gen(n,2025)[[2]]
X <- df %>% select(starts_with("X."))%>% as.matrix()

res_dir_or <- file.path(folder, "results", "oracular", "theta_opt","opt")
res_dir_est <- file.path(folder, "results", "estimated", "theta_opt")

# List full file paths
res_files_or <- list.files(path = res_dir_or, full.names = TRUE)
res_files_est <- list.files(path = res_dir_est, full.names = TRUE)

# Extract and sort by lambda number for oracular files
opt_psi_or <- as.numeric(sub(".*theta_opt_(\\d+)\\.rds", "\\1", res_files_or))
res_files_or <- res_files_or[order(opt_psi_or)]

# Extract and sort by lambda number for estimated files
opt_psi_est <- as.numeric(sub(".*theta_opt_(\\d+)\\.rds", "\\1", res_files_est))
res_files_est <- res_files_est[order(opt_psi_est)]

all_results_or <- lapply(res_files_or, readRDS)
all_psi_or <- lapply(all_results_or, make_psi)
all_psi_or_X <- lapply(all_psi_or, function(psi)psi(X))
df_all_psi_or_X <- do.call(cbind, lapply(all_psi_or_X, function(x) x))  # Assuming each entry is a vector
sigma_or <- apply(df_all_psi_or_X,2,sigma_beta)


all_results_est <- lapply(res_files_est, readRDS)
all_psi_est<- lapply(all_results_est, make_psi)
all_psi_est_X <- lapply(all_psi_est, function(psi)psi(X))
df_all_psi_est_X <- do.call(cbind, lapply(all_psi_est_X, function(x) x))  # Assuming each entry is a vector
sigma_est <- apply(df_all_psi_est_X,2,sigma_beta)


df_or <- cbind(df,sb=sigma_or)%>%as.data.frame()
df_est <- cbind(df,sb=sigma_est)%>%as.data.frame()

simu <- c("oracular","estimated")
for (type_simu in simu){
  df_std <- (if (type_simu == "oracular") df_or else df_est) %>%
  rowwise() %>%
  mutate(std_prob = sd(c_across(starts_with("sb."))))

  p1<- ggplot(df_std, aes(x = X.1, y = X.2, color = std_prob)) +
    geom_point(size = 1.2) +
    scale_color_viridis_c() +
    labs(title = "Treatment Probability Variability Over Iterations",
        color = "Std Dev") +
    theme_minimal()
  ggsave(file.path(folder,"figures",type_simu,"std_treatment_probability.pdf"),p1)


  df_mean <- (if (type_simu == "oracular") df_or else df_est) %>%
  rowwise() %>%
  mutate(mean_prob = mean(c_across(starts_with("sb."))))

  p2<- ggplot(df_mean, aes(x = X.1, y = X.2, color = mean_prob)) +
    geom_point(size = 1.2) +
    scale_color_viridis_c() +
    labs(title = "Treatment Probability Variability Over Iterations",
        color = "Mean Dev") +
    theme_minimal()

  ggsave(file.path(folder,"figures",type_simu,"mean_treatment_probability.pdf"),p2)
}

colnames(sigma_or) <- paste0("or_",1:ncol(sigma_or))
colnames(sigma_est) <- paste0("est_",1:ncol(sigma_est))

tib <- as_tibble(
    cbind(sigma_or,sigma_est)
) %>%
mutate(id=1:n()) %>% 
pivot_longer(-id, names_to="what", names_pattern="([^_]*)_.*", values_to="values") %>% 
pivot_wider(id_cols=id, names_from="what", values_from="values") %>% 
unnest(cols = c(or, est))

tib %>% 
head(n=1e5)%>%
ggplot()+
geom_point(aes(x=or,y=est),alpha=0.1, color="green")+
geom_smooth(aes(x=or,y=est, ))+
geom_abline(intercept=0,slope=1, col="red")

ggsave(file.path(folder,"figures","comparison_sigma_beta.pdf"))


########################
###### Densities #######
########################
oracular_vec <- as.vector(sigma_or)
estimated_vec <- as.vector(sigma_est)

# Combine into a data frame for plotting
densities <- data.frame(
  value = c(oracular_vec, estimated_vec),
  source = rep(c("Oracular", "Estimated"), each = length(oracular_vec))
)

# Plot density
p3 <- ggplot(densities, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Comparison of Oracular vs Estimated",
       x = "Function Output",
       y = "Density") +
  theme_minimal()

ggsave(file.path(folder,"figures","density_comparison_sigma_beta.pdf"),p3)


quantiles <- seq(0, 1, by = 0.01)

# Compute quantiles for both datasets
oracular_quantiles <- quantile(oracular_vec, probs = quantiles)
estimated_quantiles <- quantile(estimated_vec, probs = quantiles)

# Create a data frame for plotting
qq_df <- data.frame(
  quantile = rep(quantiles, 2),
  value = c(oracular_quantiles, estimated_quantiles),
  source = rep(c("Oracular", "Estimated"), each = length(quantiles))
)

# Plot ECQF (Cumulative Quantile Function)
qqdf<- ggplot(qq_df, aes(x = quantile, y = value, color = source, group = source)) +
  geom_line(size = 1.2) +
  labs(title = "ECQF: Expected Cumulative Quantile Function",
       x = "Quantile",
       y = "Function Output") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave(file.path(folder,"figures","qqdf.pdf"))

pdf(file.path(folder,"figures","qq_plot.pdf"), width = 8, height = 6)  # You can adjust the width and height

# Create the QQ plot
qqplot(oracular_vec, estimated_vec,
       main = "QQ Plot: Oracular vs Estimated",
       xlab = "Oracular Quantiles",
       ylab = "Estimated Quantiles",
       pch = 19, col = rgb(0, 0, 1, 0.5))

# Add the reference line
abline(0, 1, col = "red", lwd = 2)

# Close the PDF device (this saves the file)
dev.off()

########################
###### Densities #######
########################

source(file.path(folder,"/synthetic_data.R"))
source("src/cross-fitting.R")
library(dplyr)
library(ggplot2)
library(grf)
library(gridExtra)

df<- data_gen(n)[[2]]
X <- df %>% select(starts_with("X."))%>% as.matrix()
Jfold <-readRDS(file.path(folder,"results","data","Jfold.rds"))

folds <- CVFolds(
    N = n,
    id = NULL,
    Y = df,
    cvControl = list(V = Jfold, stratifyCV = FALSE, shuffle = TRUE)
)

folds_df <- do.call(rbind, lapply(seq_along(folds), function(v){data.frame(index=folds[[v]], fold=v)}))
s <- folds_df$fold[order(folds_df$index)]


res_dir_mu_est <- file.path(folder, "results", "data", "Mu")
res_dir_nu_est <- file.path(folder, "results", "data", "Nu")

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
ggsave(file.path(folder,"figures","CATE_std.pdf"),combined_plot)


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

combined_plot_mean <- gridExtra::grid.arrange(p21, p22, p31, p32, nrow = 2, ncol = 2)
ggsave(file.path(folder,"figures","CATE_mean.pdf"),combined_plot_mean)

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
ggsave(file.path(folder,"figures","density_comparison_CATE.pdf"),combined_plot_density)

