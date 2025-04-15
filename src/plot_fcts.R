library(gridExtra)
library(ggplot2)
library(cowplot)
library(tidyverse)


synthetic_setting_plot <- function(df_complete, delta_Y, delta_Xi){
  df_complete$sign_CATE <- as.factor(
    ifelse(
      delta_Y(
        df_complete %>% select(starts_with("X."))
      )>0,1,-1
    )
  )
  df_complete$delta_Xi <- delta_Xi(
    df_complete %>% select(starts_with("X."))
  )
    plot_Y_sign<- ggplot(df_complete, aes(x=X.1, y=X.2, color=sign_CATE))+
      geom_point(alpha = 0.5)
    p_plot<- ggplot(df_complete, aes(x=X.1, y=X.2, color=(delta_Xi)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "blue", high = "green")
    combined_plot <- grid.arrange(
      plot_Y_sign, p_plot,
      ncol = 1  
    )
  ggsave(paste0("figures/synthetic_setting.pdf"),combined_plot)
}

lambda_evol<- function(results, type_simu, beta_opt){
  results_lambda <- results[which(results$beta==beta_opt),]
  # lambda plot
  lambda_evol <- ggplot(results_lambda %>% as.data.frame(), aes(x = lambda)) +
  geom_point(aes(y = risk, color = "Risk")) +
  geom_point(aes(y = constraint, color = "Constraint")) +
  #geom_point(aes(y = policy_value, color = "Policy value")) +
  geom_point(aes(y = obj, color = "L")) +
  # Risk: with colored SE and no line
  geom_smooth(aes(y = risk, color = "Risk", fill = "Risk"), 
              method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2) +
  # Constraint: with colored SE and no line
  geom_smooth(aes(y = constraint, color = "Constraint", fill = "Constraint"), 
              method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2) +
  # Policy value: with colored SE and no line
  # geom_smooth(aes(y = policy_value, color = "Policy value", fill = "Policy value"), 
  #             method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2) +
  # L: with colored SE and no line
  geom_smooth(aes(y = obj, color = "L", fill = "L"), 
              method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2) +
  labs(
    title = expression("Evolution of optimal solution " * psi[lambda] * " with respect to " * lambda),
    subtitle = if (type_simu == "oracular") {
      expression(psi[lambda] == argmin~bgroup("{", list(L[P[0]](psi, lambda) * ": " * psi %in% Psi), "}"))
    } else {
      expression(psi[lambda] == argmin~bgroup("{", list(L[P[n]](psi, lambda) * ": " * psi %in% Psi), "}"))
    },
    x = expression(lambda),
    y = "Values"
  ) +
  # scale_color_manual(
  #   name = "Functions",
  #   values = c("Risk" = "blue", "Constraint" = "red", "L" = "gray", "Policy value" = "black"),
  #   labels = if (type_simu == "oracular") {
  #     c(
  #       "Risk" = expression(R[P[0]](psi[lambda])),
  #       "Constraint" = expression(S[P[0]](psi[lambda])),
  #       "L" = expression(L[P[0]](psi[lambda], lambda)),
  #       "Policy value" = expression(V[P[0]](pi^{"*"}))
  #     )
  #   } else {
  #     c(
  #       "Risk" = expression(R[P[n]]),
  #       "Constraint" = expression(S[P[n]]),
  #       "L" = expression(L[P[n]](psi[lambda])),
  #       "Policy value" = expression(V[P[n]](pi^{"*"}))
  #     )
  #   }
  # ) +
  # scale_fill_manual(
  #   name = "Functions",
  #   values = c("Risk" = "blue", "Constraint" = "red", "L" = "gray", "Policy value" = "black")
  # ) +
  theme_minimal()+
  guides(fill="none")
  ggsave(paste0("figures/",type_simu,"/lambda_evol.pdf"),lambda_evol)
}

sigma_beta_fixed <- function(psi_x, beta, centered){
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
    if (centered) {
        cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
    } else {
        cent <- 0
    }
    out <- c_beta * log((1 + exp(beta * psi_x)) / (1 + exp(-beta))) + cent
}

gamma_plot_funct <- function(psi_X, lambda, beta, df, centered) {
  policy<- sigma_beta_fixed(psi_X, beta, centered)
  # Initialize base plot
    p <- ggplot(
      cbind(df, treat_proba = policy),
      aes(x = X.1, y = X.2, color = treat_proba)
    ) +
    geom_point(alpha = 1) +
      scale_color_gradient2(low = "red", mid = "white", high = "blue", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish) +
      labs(title =  bquote(lambda == .(lambda) ~ "," ~ beta == .(beta))) +
      theme_minimal() +
      theme(legend.position = "right")
  return(p)
}

geom_points_fct <- function(results, idx_opt, df, type_simu, centered){
  lambda_discr <- which(results$beta==results$beta[idx_opt])[as.integer(seq(1, length(which(results$beta==results$beta[idx_opt])), length.out=10))]
  
  plots <- lapply(lambda_discr, function(x) gamma_plot_funct(results$optimal_x[[x]], results$lambda[[x]], results$beta[[x]], df, centered))
  plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))
  
  legend <- get_legend(plots[[1]])
  combined_plots <- plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
  final_plot <- plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
  # Display the final plot
  print(final_plot)
  ggsave(paste0("figures/",type_simu,"/optimal_solution_multiple_lambdas.pdf"),final_plot, width = 10, height = 6)
  
  plot_none<-gamma_plot_funct(results$optimal_x[[1]], results$lambda[[1]], results$beta[[1]], df, centered)+theme(legend.position = "none")
  plot_max <- gamma_plot_funct(results$optimal_x[[idx_opt]], results$lambda[[idx_opt]], results$beta[[idx_opt]], df, centered)+theme(legend.position = "none")

  opt_plots <- plot_grid(plot_none,plot_max, ncol=2)
  ggsave(paste0("figures/",type_simu,"/optimal_solution_optimal_lambda.pdf"),opt_plots, width = 8, height = 4)
}