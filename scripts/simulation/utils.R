library(dplyr)
library(tidyr)

# --- function to compute true parameters for each simulation setting
# df: data frame with columns "variable" and "true_value"
true_params <- function(df){
  df |>
    pivot_wider(names_from = variable, values_from = true_value) |>
    mutate(
      "delta[1]" = `shape[1]` / sqrt(1 + `shape[1]`^2),
      "delta[2]" = `shape[2]` / sqrt(1 + `shape[2]`^2),
      "mrg_scale" = sqrt(`scale[2]`^2 + `eta`^2),
      "mrg_shape" = `scale[2]` * `shape[2]` / sqrt(`scale[2]`^2 + `eta`^2 * (1 + `shape[2]`^2)),
      "mrg_delta" = `mrg_shape` / sqrt(1 + `mrg_shape`^2)
    ) |>
    mutate(
      "mean[1]" = `location[1]` + `scale[1]` * `delta[1]` * sqrt(2 / pi),
      "mean[2]" = `location[2]` + `scale[2]` * `delta[2]` * sqrt(2 / pi),
      "mrg_mean[2]" = `location[2]` + `mrg_scale` * `mrg_delta` * sqrt(2 / pi)
    ) |>
    pivot_longer(cols = 1:15, names_to = "variable", values_to = "true_value") |>
    mutate(across(where(is.numeric), round, 2))
}

# --- function to summarize posterior estimates across simulation
# df: data frame with posterior summaries for each simulation
# true_params_df: data frame with true parameter values for the simulation setting
# params2recover: vector of parameter names to recover
# setting_name: name of the simulation setting
summary_params <- function(df, true_params_df, params2recover, setting_name){
  df |>
    filter(variable %in% c(params2recover, "mu0", "mu1")) |>
    mutate(variable = factor(variable, levels = c(params2recover, "mu0", "mu1"))) |>
    arrange(sim_id, variable) |>
    group_by(sim_id) |>
    mutate(
      truth = true_params_df |> 
        filter(variable %in% c(params2recover, "mean[1]", "mrg_mean[2]")) |> 
        select(true_value) |> 
        pull(),
      bias = mean - truth,
      SE = bias^2, 
      coverage_80ci = ifelse(truth >= q10 & truth <= q90, 1, 0),
      coverage_90ci = ifelse(truth >= q5 & truth <= q95, 1, 0),
      coverage_95ci = ifelse(truth >= q2.5 & truth <= q97.5, 1, 0)
    ) |>
    ungroup() |>
    group_by(variable) |>
    summarise(
      truth = mean(truth),
      ave_mean = round(mean(mean), 2),
      avg_bias = round(mean(bias), 2),
      sd_bias = round(sd(bias), 2),
      RMSE = round(sqrt(mean(SE)),2),
      coverage_80ci = round(mean(coverage_80ci), 2),
      coverage_95ci = round(mean(coverage_95ci), 2)
    ) |>
    ungroup() |>
    mutate(setting = setting_name)
}


  
# --- function to plot s(t): annual probability of infection
api_plot_by_setting <- function(data1 = deterministic_df, 
                     data2 = spline_shat, 
                     data3 = coverage_df,
                     setting = "Setting 1"){
  ggplot() +
    # deterministic curve
    geom_line(data = data1, aes(x = time, y = s_hat),
              color = "red", size = 1.5) +
    # estimated curve
    geom_line(data = data2, aes(x = time, y = s_hat, group = sim_id),
              color = alpha("lightblue", 0.4), size = 0.7) +
    scale_x_continuous(
      breaks = seq(0, 100, by = 5)  # 5, 10, 15, 20
    ) +
    coord_cartesian(ylim = c(0, 0.105), clip = "off")  +
    labs(
      x = "",
      y = "",
      title = setting
    ) +
    theme_minimal() +                 
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid")
    ) +
    geom_text(data = data3,
              aes(x = time, y = 0.1, label = paste0(round(avg_coverage_95ci * 100,2), "%")),
              size = 2.75, angle=90,  color = "blue") + 
    annotate("text" , x = 110, y = 0.1, label = "95CI \n coverage", 
                 color = "blue", size = 2.75, fontface = "bold")
}

# --- function for analysis of classification performance 
# data_name: name of the simulated dataset
# sim_data: simulated dataset
# prob_draws: model outputs prob_noinfection[n]
eval_classification <- function(data_name, sim_data, prob_draws){
  N <- nrow(sim_data)
  # predicted probability of ever-infected = 1- prob_noinfection
  prob_inf_mean  <- 1 - colMeans(prob_draws) 
  true_status <- sim_data$Zi
  # AUC
  scores_pos <- prob_inf_mean[true_status == 1]   # predicted scores for ever-infected
  scores_neg <- prob_inf_mean[true_status == 0]   # predicted scores for never-infected
  auc <- mean(outer(scores_pos, scores_neg, ">") + 0.5 * outer(scores_pos, scores_neg, "=="))
  # Classify as ever-infected if P(ever-infected | y_n) > 0.5.
  pred_status <- as.integer(prob_inf_mean > 0.5)
  
  tp <- sum(pred_status == 1 & true_status == 1)  # true positives
  tn <- sum(pred_status == 0 & true_status == 0)  # true negatives
  fp <- sum(pred_status == 1 & true_status == 0)  # false positives
  fn <- sum(pred_status == 0 & true_status == 1)  # false negatives
  
  sensitivity <- tp / (tp + fn)   # P(predicted infected | truly infected)
  specificity <- tn / (tn + fp)   # P(predicted not infected | truly not infected)
  accuracy    <- (tp + tn) / N
  # Brier score
  brier <- mean((prob_inf_mean - true_status)^2)
  
  data.frame(
    Data = data_name,
    AUC = auc,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    Brier = brier
  )
}
