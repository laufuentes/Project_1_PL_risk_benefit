library(ggplot2)

########################
# Lambda evol #
########################

# Set up directories
res_dir_or <- file.path("MC", "results", "oracular", "group_by_lambda")
res_dir_est <- file.path("MC", "results", "estimated", "group_by_lambda")

# List full file paths
res_files_or <- list.files(path = res_dir_or, full.names = TRUE)
res_files_est <- list.files(path = res_dir_est, full.names = TRUE)

# Extract and sort by lambda number for oracular files
lambda_nums_or <- as.numeric(sub(".*lambda_(\\d+)\\.csv", "\\1", res_files_or))
res_files_or <- res_files_or[order(lambda_nums_or)]

# Extract and sort by lambda number for estimated files
lambda_nums_est <- as.numeric(sub(".*lambda_(\\d+)\\.rds", "\\1", res_files_est))
res_files_est <- res_files_est[order(lambda_nums_est)]

# Read results
all_results_or <- lapply(res_files_or, read.csv)
results_or <- do.call(rbind, all_results_or)

all_results_est <- lapply(res_files_est, readRDS)
results_est<- do.call(rbind, all_results_est)

library(tidyverse)
library(ggplot2)
library(tidyr)

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
ggsave("MC/figures/oracular/lambda_evol_or.pdf", l_evol_or)

l_evol_est<- ggplot(results_est_long, aes(x = factor(lambda), y = value, fill = variable)) + 
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +  # Separate boxplots for each variable
  labs(x = "Lambda", y = "Value", title = "Boxplots of Policy Value, Objective, Risk, and Constraint by Lambda") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Choose a color palette
ggsave("MC/figures/estimated/lambda_evol_est.pdf", l_evol_est)


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

ggsave("MC/figures/lambda_evol_both.pdf", l_evol_both)
