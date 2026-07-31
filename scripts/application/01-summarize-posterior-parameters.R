library(dplyr)

# Load model outputs
out.summary <- readRDS("data/mdl-fit-outputs-summary.rds")

# Parameters to recover
params2recover <- c(
  "location[1]", "scale[1]", "shape[1]",
  "location[2]", "scale[2]", "shape[2]",
  "eta"
)

# Check convergence of the model
out.summary |>
  filter(variable %in% c(params2recover, "mu0", "mu1" )) |>
  select(variable, rhat, ess_bulk, ess_tail) |>
  filter(rhat > 1.05) |>
  filter(ess_bulk < 400 | ess_tail < 400) 

# Summarize the posterior parameters
out.summary |>
  filter(variable %in% c(params2recover, "mu0", "mu1")) |>
  mutate(variable = factor(variable, levels = c(params2recover, "mu0", "mu1"))) |> 
  arrange(variable) |> 
  mutate(across(where(is.numeric), round, 2))

