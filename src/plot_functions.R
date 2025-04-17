library(gridExtra)
library(ggplot2)
library(cowplot)
library(tidyverse)


#' Plot Synthetic Data Setting
#'
#' Generates and saves a two-panel plot:
#' one showing the sign of the treatment effect (`delta_Mu`) and the other
#' visualizing the magnitude of selection effect (`delta_Nu`) across covariates X.1 and X.2.
#'
#' @param df_complete A data frame containing covariates prefixed with "X.".
#' @param delta_Mu A function that computes the treatment effect (mu difference) from covariates.
#' @param delta_Nu A function that computes the selection effect (nu difference) from covariates.
#'
#' @return Saves a plot to "figures/synthetic_setting.pdf".
#' @export
synthetic_setting_plot <- function(df_complete, delta_Mu, delta_Nu) {
  `%>%` <- magrittr::`%>%`
  df_complete$sign_delta_Mu <- as.factor(
    sign(
      delta_Mu(df_complete %>% dplyr::select(dplyr::starts_with("X.")))
    )
  )
  df_complete$delta_Nu <- delta_Nu(
    df_complete %>% dplyr::select(dplyr::starts_with("X."))
  )
  plot_Y_sign <- ggplot2::ggplot(df_complete, ggplot2::aes(x = df_complete$X.1, y = df_complete$X.2, color = df_complete$sign_delta_Mu)) +
    ggplot2::geom_point(alpha = 0.5)
  p_plot <- ggplot2::ggplot(df_complete, ggplot2::aes(x = df_complete$X.1, y = df_complete$X.2, color = (df_complete$delta_Nu))) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_gradient(low = "blue", high = "green")
  combined_plot <- gridExtra::grid.arrange(
    plot_Y_sign, p_plot,
    ncol = 1
  )
  ggplot2::ggsave(file.path("figures", "synthetic_setting.pdf"), combined_plot)
}

#' Plot Evolution of Objective Terms Across Lambda Values
#'
#' Visualizes how risk, constraint, and overall objective evolve with respect to different lambda values.
#' Includes smooth loess trends and confidence intervals.
#'
#' @param results A data frame or tibble containing `lambda`, `risk`, `constraint`, `obj`, and `beta` columns.
#' @param type_simu A string indicating the simulation type (e.g., "oracular" or "empirical").
#' @param beta_opt A numeric value for the beta whose lambda evolution is being plotted.
#'
#' @return Saves a plot to "figures/<type_simu>/lambda_evol.pdf".
#' @export
lambda_evol <- function(results, type_simu, beta_opt) {
  results_lambda <- results[which(results$beta == beta_opt), ]
  res <- as.data.frame(results_lambda)
  # lambda plot
  lambda_evol_plot <- ggplot2::ggplot(res, ggplot2::aes(x = res$lambda)) +
    ggplot2::geom_point(ggplot2::aes(y = res$risk, color = "Risk")) +
    ggplot2::geom_point(ggplot2::aes(y = res$constraint, color = "Constraint")) +
    #geom_point(aes(y = policy_value, color = "Policy value")) +
    ggplot2::geom_point(ggplot2::aes(y = res$obj, color = "L")) +
    # Risk: with colored SE and no line
    ggplot2::geom_smooth(ggplot2::aes(y = res$risk, color = "Risk", fill = "Risk"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    # Constraint: with colored SE and no line
    ggplot2::geom_smooth(ggplot2::aes(y = res$constraint, color = "Constraint", fill = "Constraint"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    ggplot2::geom_smooth(ggplot2::aes(y = res$obj, color = "L", fill = "L"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    ggplot2::labs(
      title = expression("Evolution of optimal solution " * psi[lambda] * " with respect to " * lambda),
      subtitle = if (type_simu == "oracular") {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[0]](psi, lambda) * ": " * psi %in% Psi), "}"))
      } else {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[n]](psi, lambda) * ": " * psi %in% Psi), "}"))
      },
      x = expression(lambda),
      y = "Values"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(fill = "none")
  
  ggplot2::ggsave(file.path("figures", type_simu, "lambda_evol.pdf"), lambda_evol_plot)
}

#' Visualize Treatment Assignment Probability
#'
#' Plots the smoothed treatment assignment probability over covariates X.1 and X.2.
#'
#' @param psi_X A numeric vector of policy scores.
#' @param lambda The lambda value used (for display purposes).
#' @param beta The smoothing parameter used in the policy.
#' @param df A data frame with at least columns X.1 and X.2.
#' @param centered Logical indicating whether probabilities are centered.
#'
#' @return A ggplot object.
#' @export
gamma_plot_funct <- function(psi_X, lambda, beta, df, centered) {
  policy <- sigma_beta(psi_X, beta, centered)
  # Initialize base plot
  p <- ggplot2::ggplot(
    cbind(df, treat_proba = policy),
    ggplot2::aes(x = df$X.1, y = df$X.2, color = df$treat_proba)
  ) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::scale_color_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0.5,
      limits = c(0, 1),
      oob = scales::squish
    ) +
    ggplot2::labs(title = bquote(lambda == .(lambda) ~ "," ~ beta == .(beta))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  return(p)
}

#' Plot Optimal Policy Across Discrete Lambda Values
#'
#' Generates multiple plots of treatment assignment probability for a range of lambda values
#' and combines them with a shared legend. Also saves plots for the initial and optimal lambdas.
#'
#' @param results A list or data frame containing `optimal_x`, `lambda`, `beta`, etc.
#' @param idx_opt The index of the optimal lambda in the results.
#' @param df A data frame with covariates, including X.1 and X.2.
#' @param type_simu A string used to name the output folder (e.g., "oracular", "estimated").
#' @param centered Logical; if TRUE, uses centered sigma transformation.
#'
#' @return Saves plots in "figures/<type_simu>/".
#' @export
geom_points_fct <- function(results, idx_opt, df, type_simu, centered){
  lambda_discr <- which(results$beta==results$beta[idx_opt])[as.integer(seq(1, length(which(results$beta==results$beta[idx_opt])), length.out=10))]
  
  plots <- lapply(lambda_discr, function(x) gamma_plot_funct(results$optimal_x[[x]], results$lambda[[x]], results$beta[[x]], df, centered))
  plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
  
  legend <- cowplot::get_legend(plots[[1]])
  combined_plots <- cowplot::plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
  final_plot <- cowplot::plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
  # Display the final plot
  print(final_plot)
  ggplot2::ggsave(file.path("figures",type_simu,"optimal_solution_multiple_lambdas.pdf"),final_plot, width = 10, height = 6)
  
  plot_none<-gamma_plot_funct(results$optimal_x[[1]], results$lambda[[1]], results$beta[[1]], df, centered)+ggplot2::theme(legend.position = "none")
  plot_max <- gamma_plot_funct(results$optimal_x[[idx_opt]], results$lambda[[idx_opt]], results$beta[[idx_opt]], df, centered)+ggplot2::theme(legend.position = "none")

  opt_plots <- cowplot::plot_grid(plot_none,plot_max, ncol=2)
  ggplot2::ggsave(file.path("figures",type_simu,"optimal_solution_optimal_lambda.pdf"),opt_plots, width = 8, height = 4)
}