# Simulate Lao malaria serology dataset
# log(PvAMA1) ~ Skew Normal mixture distribution

# Parameters ----------------------------------------------------------------
# times: number of time points
# birth_times: vector of birth times of n individuals
# prob_infection_fun: a smooth function of infection rates over time
# NOTE: if include_randeff = TRUE, then random effects are included
# only for the infected individuals.
# waning: logical, whether to include waning effects (default FALSE)
# eta: standard deviation of random effects (default 1)
# rho: antibody waning rate (default 0)
# xi1, omega1, alpha1: parameters of skew-normal distribution for infected individuals
# xi0, omega0, alpha0: parameters of skew-normal distribution for never infected individuals


simulate_serodata <- function(times, birth_times,
                              prob_infection_fun,
                              waning = FALSE,
                              eta, rho = 0, 
                              xi0, omega0, alpha0,
                              xi1, omega1, alpha1) {
  # base data frame with infection status
  base_dat <- data.frame(birth_times = birth_times) |>
    mutate(
      age = max(times) - birth_times
    ) |>
    rowwise() |>
    mutate(
      # calculate probability of never being infected by the age
      prob_never_infected = prod(1 - prob_infection_fun(max(times) - 0:age)),
      # infection status (Z_i=1: infected, Z_i=0: not infected)
      Zi = rbinom(1, 1, 1 - prob_never_infected)
    ) |>
    # label infection status
    mutate(
      infection_status = factor(Zi,
        levels = c(0, 1),
        labels = c("Not infected", "Infected")
      )
    ) |>
    ungroup()

  # determine infected_mean based on infection status and random effects
  u_i <- rnorm(nrow(base_dat), mean = 0, sd = 1)
  if (!waning) {
    # no waning
    dat <- base_dat |>
      mutate(xi_1i = ifelse(infection_status == "Infected", xi1 + eta * u_i, NA))
  } else {
    # with random effects and waning
    dat <- base_dat |>
      rowwise() |>
      mutate(
        time_infected = ifelse(infection_status == "Infected",
          sample(1:age, 1), NA
        ),
        waining_time = age - time_infected,
        xi_1i = ifelse(infection_status == "Infected", 
                       xi1 + eta * u_i + rho * waining_time, NA)
      ) |>
      ungroup() |>
      select(-time_infected, -waining_time)
  } 

  # log_PvAMA1
  df <- dat |>
    rowwise() |>
    mutate(
      log_PvAMA1 = ifelse(
        infection_status == "Infected",
        sn::rsn(1, xi = xi_1i, omega = omega1, alpha = alpha1),
        sn::rsn(1, xi = xi0, omega = omega0, alpha = alpha0)
      )
    ) |>
    ungroup() |>
    select(-xi_1i)

  out <- df |> mutate(PvAMA1 = exp(log_PvAMA1))

  return(out)
}
