library(dplyr)
source("scripts/simulation/utils.R")
# Load model fit results from setting 1
df_results_setting1 <- readRDS("data/mdl-fit-simdata-setting1-outputs-summary.rds")
df_results_setting2 <- readRDS("data/mdl-fit-simdata-setting2-outputs-summary.rds")
df_results_setting3 <- readRDS("data/mdl-fit-simdata-setting3-outputs-summary.rds")

# Define parameters to recover
params2recover <- c(
  "location[1]", "scale[1]", "shape[1]",
  "location[2]", "scale[2]", "shape[2]",
  "eta"      
)

# true parameter values
df1 <- data.frame(
  variable = params2recover,
  true_value = c(3.5, 0.5, -3.5, 5.5, 1, 3.5, 0)
)
true_params_setting1 <- true_params(df1)

df2 <- data.frame(
  variable = params2recover,
  true_value = c(4.2, 1, -3.5, 5, 1.5, 3.5, 0)
) 
true_params_setting2 <- true_params(df2)

df3 <- data.frame(
  variable = params2recover,
  true_value = c(4.2, 1, -3.5, 5, 1.5, 3.5, 1)
)
true_params_setting3 <- true_params(df3)
  
  
# summarize posterior estimates across simulations
tbl_summary_params_setting1 <- summary_params(
  df_results_setting1, true_params_setting1, params2recover, "setting1"
  )

tbl_summary_params_setting2 <- summary_params(
  df_results_setting2, true_params_setting2, params2recover, "setting2"
)

tbl_summary_params_setting3 <- summary_params(
  df_results_setting3, true_params_setting3, params2recover, "setting3"
)

# Combine summaries from all settings into a single table
summary_all <- tbl_summary_params_setting1 |>
  bind_rows(tbl_summary_params_setting2) |>
  bind_rows(tbl_summary_params_setting3) |>
  select(setting, variable, truth, ave_mean, avg_bias, sd_bias, RMSE, coverage_80ci, coverage_95ci)
