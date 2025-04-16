# load_results.R
library(tidyverse)
centered <- FALSE
df <- read.csv(file.path("results","estimated","df.csv"), stringsAsFactors=FALSE)

 # Directory where individual result files are saved
 res_dir <- file.path("results","estimated","individual_results")
 
 # List all result files
 res_files <- list.files(res_dir, pattern = "^res_\\d+\\.rds$", full.names = TRUE)
 
 # Load all results into a list
 all_results <- lapply(res_files, readRDS)
 
 # Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
 # This step depends on what `res` actually is in your process
 # If each `res` is a named list or a data.frame row:
 results<- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    lambda = x[[1]],
    beta = x[[2]],
    optimal_x = I(x[[3]]), # Take only the first element of optimal_x
    risk = x[[4]],
    constraint = x[[5]],
    obj = x[[6]],
    policy_value = x[[7]]
  )
}))
 
 # Save combined results
saveRDS(results, file = file.path("results","estimated","results.rds"))
write.csv(results%>%select(-optimal_x), file = file.path("results","estimated","results.csv"), row.names = FALSE)

result_T <- results %>%
  group_by(lambda, beta) %>%
  summarise(
    obj = mean(obj),
    risk = mean(risk),
    constraint = mean(constraint),
    policy_value = mean(policy_value),
    .groups = "drop"
  )

write.csv(
  result_T,
  file.path("results","estimated","result_T.csv")
)

saveRDS(result_T, file = file.path("results","estimated","result_T.rds"))


idx_opt_obj <- which(
  result_T$obj == max(result_T$obj[result_T$constraint <= 0])
)
idx_opt <- idx_opt_obj

source(file.path("src","plot_functions.R"))
lambda_evol(
    result_T, 
    "estimated", 
    results$beta[[idx_opt]]
)

source(file.path("src","algorithm_architecture.R"))

theta_lambda0 <- Final_policy(0, results$beta[[idx_opt]],X, s,  Delta_mu_nj_folds, Delta_nu_nj_folds)
psi0<- make_psi(theta_lambda0)
psi_x0 <- psi0(X)

theta_opt <- Final_policy(lambda, results$lambda[[idx_opt]],X, s, Delta_mu_nj_folds, Delta_nu_nj_folds)
optimal_psi <- make_psi(theta_opt)
optimal_x <- optimal_psi(X)

p0 <- gamma_plot_funct(psi_x0, 0, 0.05, df, centered)+ theme(legend.position = "none")
p <- gamma_plot_funct(optimal_x, lambda, 0.05, df, centered)+ theme(legend.position = "none")
opt_plots <- plot_grid(p0,p, ncol=2)
ggsave(file.path("figures","estimated","Optimal_policy.pdf"),opt_plots, width = 8, height = 4)

# sb_est <- sigma_beta(optimal_psi, X, 0.05, centered)
# sb_or <- sigma_beta()

# tibble(
#   est=sb_est, 
#   ora=sb_or
# )%>%
# ggplot()+
# geom_point(aes(x=est,y=ora))+
# geom_abline(intercept=0,slope=1, col="red")

# ggsave("images/estimation/toto.pdf")