source("src/optim_functions.R")
source("src/utils.R")
source("src/synthetic_data.R")

library(dplyr)
library(tidyverse)
library(ggplot2)

########################
# Psi opt #
########################

df <- data_gen(1e4,2025)[[2]]
X <- df %>% select(starts_with("X."))%>% as.matrix()

res_dir_or <- file.path("MC", "results", "oracular", "theta_opt","opt")
res_dir_est <- file.path("MC", "results", "estimated", "theta_opt")

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

type_simu <- "oracular"

df_std <- (if (type_simu == "oracular") df_or else df_est) %>%
  rowwise() %>%
  mutate(std_prob = sd(c_across(starts_with("sb."))))

p1<- ggplot(df_std, aes(x = X.1, y = X.2, color = std_prob)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c() +
  labs(title = "Treatment Probability Variability Over Iterations",
       color = "Std Dev") +
  theme_minimal()

ggsave(file.path("MC","figures",type_simu,"std_treatment_probability.pdf"),p1)


df_mean <- (if (type_simu == "oracular") df_or else df_est) %>%
  rowwise() %>%
  mutate(mean_prob = mean(c_across(starts_with("sb."))))

p2<- ggplot(df_mean, aes(x = X.1, y = X.2, color = mean_prob)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c() +
  labs(title = "Treatment Probability Variability Over Iterations",
       color = "Std Dev") +
  theme_minimal()

ggsave(file.path("MC","figures",type_simu,"mean_treatment_probability.pdf"),p2)


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

ggsave(file.path("MC","figures","density_comparison.pdf"),p3)


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
ggplot(qq_df, aes(x = quantile, y = value, color = source, group = source)) +
  geom_line(size = 1.2) +
  labs(title = "ECQF: Expected Cumulative Quantile Function",
       x = "Quantile",
       y = "Function Output") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("try.pdf")

pdf("qq_plot.pdf", width = 8, height = 6)  # You can adjust the width and height

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
# CATEs#
########################


