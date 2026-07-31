library(dplyr)

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

# Identify simulations with convergence issues
nonconv_sim_id_setting1 <- df_results_setting1 |>
  filter(variable %in% c(params2recover, "mu0", "mu1" )) |>
  select(sim_id, variable, rhat, ess_bulk, ess_tail) |>
  filter(rhat > 1.05) |>
  filter(ess_bulk < 400 | ess_tail < 400) |>
  pull(sim_id) |>
  unique()

nonconv_sim_id_setting2 <- df_results_setting2 |>
  filter(variable %in% c(params2recover, "mu0", "mu1" )) |>
  select(sim_id, variable, rhat, ess_bulk, ess_tail) |>
  filter(rhat > 1.05) |>
  filter(ess_bulk < 400 | ess_tail < 400) |>
  pull(sim_id) |>
  unique()

nonconv_sim_id_setting3 <- df_results_setting3 |>
  filter(variable %in% c(params2recover, "mu0", "mu1" )) |>
  select(sim_id, variable, rhat, ess_bulk, ess_tail) |>
  filter(rhat > 1.05) |>
  filter(ess_bulk < 400 | ess_tail < 400) |>
  pull(sim_id) |>
  unique()

