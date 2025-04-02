setwd("~/Project_1_PL_risk_benefit")

set.seed(2025)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(gganimate)

source("src/synthetic_data.R")

synthetic_setting_plot <- function(df_complete){
  ## Create a the variable for plotting CATE-Y
  df_complete$sign_CATE <- as.factor(
    ifelse(
      delta_Y(
        df_complete %>% select(starts_with("X"))
      )>0,1,-1
    )
  )

  ## Create the CATE-Z variable
  df_complete$delta_Xi <- delta_Xi(
    df_complete %>% select(starts_with("X"))
  )
  
  # Plot
  if(option[1]=="IVF"){
    plot_Y_sign<- ggplot(df_complete, aes(x=X.height, y=X.bmi, color=sign_CATE))+
      geom_point(alpha = 0.5)
    p0_plot<- ggplot(df_complete, aes(x=X.height, y=X.bmi, color=(p0)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "yellow", high = "red", limits = c(0, 1))
    p1_plot<- ggplot(df_complete, aes(x=X.height, y=X.bmi, color=(p1)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "yellow", high = "red", limits = c(0, 1))
    
    p_plot<- ggplot(df_complete, aes(x=X.height, y=X.bmi, color=(delta_Xi)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "blue", high = "green", limits = c(0.1, 0.6))
    
    if (option[2] == "2") {
      p_plot <- p_plot +
        geom_rect(aes(xmin = min(X.height), xmax = 150, ymin = min(X.bmi), ymax = 20), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
        geom_rect(aes(xmin = min(X.height), xmax = 155, ymin = min(X.bmi), ymax = 25), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
    }
    combined_plot <- grid.arrange(
      plot_Y_sign, p0_plot, p1_plot, p_plot,
      ncol = 1  
    )
  }else{
    plot_Y_sign<- ggplot(df_complete, aes(x=X.1, y=X.2, color=sign_CATE))+
      geom_point(alpha = 0.5)
    p0_plot<- ggplot(df_complete, aes(x=X.1, y=X.2, color=(Xi_0)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "yellow", high = "red", limits = c(0,1))
    p1_plot<- ggplot(df_complete, aes(x=X.1, y=X.2, color=(Xi_1)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "yellow", high = "red", limits = c(0, 1))
    
    p_plot<- ggplot(df_complete, aes(x=X.1, y=X.2, color=(delta_Xi)))+
      geom_point(alpha = 0.5)+
      scale_color_gradient(low = "blue", high = "green")
    combined_plot <- grid.arrange(
      plot_Y_sign, p0_plot, p1_plot, p_plot,
      ncol = 1  
    )
  }
  ggsave(paste0("images/synthetic_setting_",option[1],"_",option[2],".pdf"),combined_plot)
}

lambda_evol<- function(results, type_simu, option, lambda_opt, beta_opt){
  results_beta <- results[which(results$lambda==lambda_opt),]
  results_lambda <- results[which(results$beta==beta_opt),]

  # lambda plot
  # lambda plot
  lambda_evol <- ggplot(results_lambda %>% as.data.frame(), aes(x = lambda)) +
    geom_point(aes(y = risk, color = "Risk")) +
    geom_point(aes(y = constraint, color = "Constraint")) +
    geom_point(aes(y = policy_value, color = "Policy value")) +
    geom_point(aes(y = obj, color = "L")) +
    geom_smooth(aes(y = risk, color = "Risk"), method = "loess", se = FALSE, span = 0.5) +
    geom_smooth(aes(y = constraint, color = "Constraint"), method = "loess", se = FALSE, span = 0.5) +
    geom_smooth(aes(y = policy_value, color = "Policy value"), method = "loess", se = FALSE, span = 0.5) +
    geom_smooth(aes(y = obj, color = "L"), method = "loess", se = FALSE, span = 0.5) +
    labs(
      title= expression("Evolution of optimal solution " * psi[lambda]*" with respect to " * lambda), # Added expression here
      subtitle = if (type_simu == "oracular") {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[0]](psi, lambda) * ": " * psi %in% Psi), "}"))
      } else {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[n]](psi, lambda) * ": " * psi %in% Psi), "}"))
      },
      x = expression(lambda), # added expression here.
      y = "Values"
    ) +
    scale_color_manual(
    name = "Functions",
    values = c("Risk" = "blue", "Constraint" = "red", "L" = "gray", "Policy value" = "black"),
    labels = if (type_simu == "oracular") {
      c(
        "Risk" = expression(R[P[0]](psi[lambda])),
        "Constraint" = expression(S[P[0]](psi[lambda])),
        "L" = expression(L[P[0]](psi[lambda],lambda)),
        "Policy value" = expression(V[P[0]](pi^{"*"}))
      ) # Closing parenthesis for c()
    } else {
      c(
        "Risk" = expression(R[P[n]]),
        "Constraint" = expression(S[P[n]]),
        "L" = expression(L[P[n]](psi[lambda])),
        "Policy value" = expression(V[P[n]](pi^{"*"}))
      ) # Closing parenthesis for c()
    }
  )+theme_minimal()
  ggsave(paste0("images/",type_simu,"/lambda_evol_",option[1],"_",option[2],".pdf"),lambda_evol)

  # beta plot
  # beta_evol<- ggplot(results_beta, aes(x = beta)) +
  #   geom_point(aes(y = risk, color = "Risk")) +
  #   geom_point(aes(y = constraint, color = "Constraint")) +
  #   geom_point(aes(y=obj, color="L"))+
  #   geom_point(aes(y=policy_value,color="Policy value"))+
  #   geom_line(aes(y = risk, color = "Risk")) +
  #   geom_line(aes(y = constraint, color = "Constraint")) +
  #   geom_line(aes(y=policy_value,color="Policy value"))+
  #   geom_line(aes(y=obj, color="L"))+
  #   labs(title = "Risk and Constraint vs Beta", x = "Beta", y = "Value") +
  #   scale_color_manual(name = "Functions", values = c("Risk" = "blue", "Constraint" = "red", "L"="gray", "Policy value"="black")) +
  #   theme_minimal()
  
  # ggsave(paste0("images/",type_simu,"/beta_evol_",option[1],"_",option[2],".pdf"),beta_evol)
}

policy_value_plot<-  function(results,type_simu, lambda_opt, beta_opt, option){
  results_beta <- results[which(results$lambda==lambda_opt),]
  results_lambda <- results[which(results$beta==beta_opt),]
  # lambda plot
  lambda_evol<- ggplot(results_lambda %>% as.data.frame(), aes(x = lambda)) +
    geom_point(aes(y=policy_value, color="Policy value"))+
    geom_point(aes(y=obj, color="L"))+
    labs(title = "Policy value vs Lambda", x = "Lambda", y = "Policy Value") +
    scale_color_manual(name = "Functions", values = c("Policy value" = "black","L" = "gray")) +
    theme_minimal()
  
  ggsave(paste0("images/",type_simu,"/lambda_evol_",option[1],"_",option[2],".pdf"),lambda_evol)
  
  # beta plot
  beta_evol<- ggplot(results_beta, aes(x = beta)) +
    geom_point(aes(y = risk, color = "Risk")) +
    geom_point(aes(y = constraint, color = "Constraint")) +
    geom_point(aes(y=obj, color="L"))+
    labs(title = "Risk and Constraint vs Beta", x = "Beta", y = "Value") +
    scale_color_manual(name = "Functions", values = c("Risk" = "blue", "Constraint" = "red", "L"="gray")) +
    theme_minimal()
  
  ggsave(paste0("images/",type_simu,"/beta_evol_",option[1],"_",option[2],".pdf"),beta_evol)
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

gamma_plot_funct <- function(psi_X, lambda, beta, df, option, centered) {
  policy<- sigma_beta_fixed(psi_X, beta, centered)
  # Initialize base plot
  if(option[1]=="IVF"){
    p <- ggplot(
      cbind(df, treat_proba = policy_FW),
      aes(x = X.height, y = X.bmi, color = treat_proba)
    ) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish) +
      labs(title =  bquote(lambda == .(lambda) ~ "," ~ beta == .(beta))) +
      theme_minimal() +
      theme(legend.position = "right")
    
    # Add rectangles **only if option == 1**
    if (option[2] == 2) {
      p <- p +
        geom_rect(aes(xmin = min(X.height), xmax = 150, ymin = min(X.bmi), ymax = 20), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
        geom_rect(aes(xmin = min(X.height), xmax = 155, ymin = min(X.bmi), ymax = 25), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
    }
  }else{
    p <- ggplot(
      cbind(df, treat_proba = policy),
      aes(x = X.1, y = X.2, color = treat_proba)
    ) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish) +
      labs(title =  bquote(lambda == .(lambda) ~ "," ~ beta == .(beta))) +
      theme_minimal() +
      theme(legend.position = "right")
  }
  
  return(p)
}

geom_points_fct <- function(results, idx_opt, df, type_simu, option, centered=TRUE){
  lambda_discr <- which(results$beta==results$beta[idx_opt])[as.integer(seq(1, length(which(results$beta==results$beta[idx_opt])), length.out=10))]
  
  plots <- lapply(lambda_discr, function(x) gamma_plot_funct(results$optimal_x[[x]], results$lambda[[x]], results$beta[[x]], df,option, centered))
  plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))
  
  legend <- get_legend(plots[[1]])
  combined_plots <- plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
  final_plot <- plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
  # Display the final plot
  print(final_plot)
  ggsave(paste0("images/",type_simu,"/geom_points_lambda_",option[1],"_",option[2],".pdf"),final_plot, width = 16, height = 10)
  
  plot_none<-gamma_plot_funct(results$optimal_x[[1]], results$lambda[[1]], results$beta[[1]], df,option, centered)
  plot_min <- gamma_plot_funct(results$optimal_x[[idx_opt-1]], results$lambda[[idx_opt-1]], results$beta[[idx_opt-1]], df,option, centered)
  plot_max <- gamma_plot_funct(results$optimal_x[[idx_opt]], results$lambda[[idx_opt]], results$beta[[idx_opt]], df,option, centered)

  opt_plots <- plot_grid(plot_none, plot_min,plot_max, ncol=3)
  ggsave(paste0("images/",type_simu,"/geom_points_lambda_opt_",option[1],"_",option[2],".pdf"),opt_plots, width = 16, height = 8)

}

gamma_lambda_plot <- function(results, df, type_simu, option, centered=TRUE){
  # Select the row with the maximum policy_value for each lambda
  lambda_results <- results %>% 
    select(-optimal_x) %>% 
    group_by(lambda) %>%
    filter(constraint<0) %>% 
    slice_max(policy_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>% 
    as.data.frame()
  
  lambda_plot <- ggplot(lambda_results, aes(lambda))+
    geom_point(aes(y=beta))+
    geom_line(aes(y=beta))+
    theme_minimal()
  
  ggsave(paste0("images/",type_simu,"/lambda_plot_",option[1],"_",option[2],".pdf"),lambda_plot)
  
  beta_res <- results %>%
    group_by(beta) %>%
    filter(constraint<0) %>% 
    slice_max(policy_value, n = 1, with_ties = FALSE) %>%
    ungroup()
    
  beta_results <- beta_res %>% 
    select(-optimal_x) %>% 
    as.data.frame()
  
  beta_plot <-ggplot(beta_results, aes(beta, y=lambda))+
    geom_point()+
    geom_line()+
    theme_minimal()
  
  ggsave(paste0("images/",type_simu,"/beta_plot_",option[1],"_",option[2],".pdf"),beta_plot)
  
  beta_discr <- as.integer(seq(1, dim(beta_res)[1], length.out=10))
  plots <- lapply(beta_discr, function(x) gamma_plot_funct(beta_res$optimal_x[[x]], beta_res$lambda[[x]], beta_res$beta[[x]], df,option, centered))
  plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))
  
  legend <- get_legend(plots[[1]])
  combined_plots <- plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
  final_plot_b <- plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
  print(final_plot_b)
  
  ggsave(paste0("images/",type_simu,"/geom_beta_plot_",option[1],"_",option[2],".pdf"),final_plot_b)
  
}

animation_plot<- function(results, lambda_values, df, type_simu, option, centered){
  # For lambda 
  if(option[1]=="IVF"){
    df_covariates <- df %>% 
      select("X.bmi","X.height") %>% 
      mutate(id=1:n)
    
    df_policy <- as.data.frame(sapply(results$optimal_x, function(x){sigma_beta(x, beta, centered)}))
    colnames(df_policy) <- paste0("lambda_",lambda_values)  # Adjust step size if needed
    df_policy$id <- 1:n
    
    df_long <- df_policy %>%
      pivot_longer(cols = starts_with("lambda_"),
                   names_to = "lambda",
                   values_to = "policy") %>%
      mutate(lambda = as.numeric(sub("lambda_", "", lambda)))
    
    df_final <- df_long %>%
      left_join(df_covariates, by = "id")
    
    
    p <- ggplot(df_final, aes(x = X.height, y = X.bmi, color = policy)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish)+
      labs(title = "Treatment Probability for λ = {closest_state}",
           x = "Height", y = "BMI", color = "Treatment Probability") +
      transition_states(lambda) +
      ease_aes('linear')
    
    if (option[2] == 2) {
      p <- p +
        geom_rect(aes(xmin = min(X.height), xmax = 150, ymin = min(X.bmi), ymax = 20), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
        geom_rect(aes(xmin = min(X.height), xmax = 155, ymin = min(X.bmi), ymax = 25), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
    }
  }
  else{
    df_covariates <- df %>% 
      select("X.1","X.2") %>% 
      mutate(id=1:n)
    
    df_policy <- as.data.frame(sapply(results$optimal_x, function(x){sigma_beta(x, beta, centered)}))
    colnames(df_policy) <- paste0("lambda_",lambda_values)  # Adjust step size if needed
    df_policy$id <- 1:n
    
    df_long <- df_policy %>%
      pivot_longer(cols = starts_with("lambda_"),
                   names_to = "lambda",
                   values_to = "policy") %>%
      mutate(lambda = as.numeric(sub("lambda_", "", lambda)))
    
    df_final <- df_long %>%
      left_join(df_covariates, by = "id")
    
    
    p <- ggplot(df_final, aes(x = X.1, y = X.2, color = policy)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish)+
      labs(title = "Treatment Probability for λ = {closest_state}",
           x = "X1", y = "X2", color = "Treatment Probability") +
      transition_states(lambda) +
      ease_aes('linear')
    
  }
  animate(p, fps = 20, duration = 10, width = 800, height = 600)
  anim_save(paste0("images/",type_simu,"/policy_animation_",option[1],"_",option[2],".gif"))
  
}
