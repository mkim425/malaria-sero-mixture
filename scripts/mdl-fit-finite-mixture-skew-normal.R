library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)
library(sn) #package for skew-normal dist.
library(purrr)
library(dplyr)
library(shinystan)
options(mc.cores = 4)


data <- readRDS("...") # add appropriate file name
N <- nrow(data)

time_points <- 1:100
num_times <- length(time_points)

# Create B-spline basis matrix
B <- bs(time_points, 
        df = 6, 
        intercept = TRUE)
num_basis <- ncol(B) # (internal knots + degree = 6)


# prepare model
mdl_SN_Pspline <- cmdstan_model(
  stan_file = "finite-mixture-skew-normal.stan")


# initial values function for all chains to avoid pathological starts
init_fun <- function(chain_id = 1) {
  list(
    location = c(3.5, 5.8),   # near prior means; ensures location[1] < location[2]
    scale    = c(0.8, 1.2),   # near half-t(3, 0, 1.5) prior median (~1.15)
    shape    = c(0, 0),       # start with symmetric (no skew) for both components
    eta      = 0.2,           # reasonable positive start for half-Cauchy(0, 1); satisfies lower bound 0
    c1       = -3,            # logit(s) = -3 => s ~ 0.047 (low baseline infection risk)
    c2       = -3,
    z_k      = rep(0, num_basis - 2),  # zero innovations -> flat linear spline
    tau      = 0.2            # moderate smoothness for s(t) spline
  )
}

# Fit the model with transformed parameters of B-spline with penalized priors
fit_mdl_Pspline <- mdl_SN_Pspline$sample(
  data = list(N = N,
              y = log(data$PvAMA1),
              K = 2,
              a = data$age,
              num_basis = num_basis,
              num_times = num_times,
              B = B),
  seed = 12345,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 3000,
  adapt_delta = 0.999,
  max_treedepth = 15,
  init = init_fun, 
  sig_figs = 6
)

# Diagnostics: divergent transitions and max treedepth saturation
diag_sum <- fit_mdl_Pspline$diagnostic_summary()
# Reported quantities:
#   num_divergent    : divergent transitions per chain (target: 0)
#   num_max_treedepth: transitions hitting max_treedepth (target: 0)
#   ebfmi            : E-BFMI per chain (target: > 0.3; low values suggest funnel)
print(diag_sum)


# Output: parameter summary and posterior draws
out1 <- fit_mdl_Pspline$summary()
out2 <- fit_mdl_Pspline$draws() |>
  summarise_draws(
    ~quantile(
    .x, 
    probs = c(0.01, 0.025, 0.10, 0.25, 0.50, 0.75, 0.9, 0.975, 0.99)
    )
  ) |>
  rename(
    q1 = `1%`,
    q2.5 = `2.5%`,
    q10 = `10%`,
    q25 = `25%`,
    q50 = `50%`,
    q75 = `75%`,
    q90 = `90%`,
    q97.5 = `97.5%`,
    q99 = `99%`
  )

# output: parameter summary
summary_df <- out1 |>
  left_join(out2, by = "variable") |>
  mutate(sim_id = file_name)

# save parameter summary
saveRDS(summary_df, file = "../output/mdl-fit-summary")

# save posterior draws
draws <- fit_mdl_Pspline$draws()
saveRDS(draws, "../output/mdl-fit-draws")
