# generate simulated data: 1000 datasets, each containing 1,000 individuals for each setting

library(dplyr)
library(sn) # for skew normal distribution

source("./scripts/simulation/annual_prob_infection.R")
#==> function annual_prob_infection()
source("./scripts/simulation/simulate_serodata.R")
#==> function: simulate_serodata(n, times, birth_times,
#                               prob_infection_fun,
#                               waning = FALSE,
#                               eta = 1, 
#                               xi0, omega0, alpha0, xi1, omega1, alpha1)

### load training data ------------------------------------------------------
training_data <- read.csv("./data/training_serology_sample.csv")

# annual probability of infection function
s <- annual_prob_infection(0.00055, 0.055, 2, 0.1, 70)
prob_infection_fun <- s

# generate birth times based on age distribution in the training data
sample_size <- 1000
age_distribution <- table(training_data$AGE) / nrow(training_data) # age distribution in the training data
times <- seq(1, 100, by = 1) # simulate yearly time steps across a 100 year period

# simulate 1000 datasets for each setting
for (i in 1:1000) {
  # sample ages based on the age distribution in the training data
  set.seed(i)
  sampled_ages <- sample(
    names(age_distribution),
    size = sample_size,
    replace = TRUE,
    prob = age_distribution
  ) |>
    as.numeric()
  # calculate birth times
  birth_times <- max(times) - sampled_ages

  # setting 1:
  set.seed(1000+i) # set seed for reproducibility
  df1_sn <- simulate_serodata(
    times, birth_times, prob_infection_fun,
    waning = FALSE,
    eta = 0, 
    xi0 = 3.5, omega0 = 0.5, alpha0 = -3.5,
    xi1 = 5.5, omega1 = 1, alpha1 = 3.5
  )

  # setting 2:
  set.seed(2000+i) # set seed for reproducibility
  df2_sn <- simulate_serodata(
    times, birth_times, prob_infection_fun,
    waning = FALSE,
    eta = 0, 
    xi0 = 4.2, omega0 = 1, alpha0 = -3.5,
    xi1 = 5, omega1 = 1.5, alpha1 = 3.5
  )

  # setting 3: random effect ~ Normal(0, 1)
  set.seed(3000+i) # set seed for reproducibility
  df3_sn <- simulate_serodata(
    times, birth_times, prob_infection_fun,
    waning = FALSE,
    eta = 1, 
    xi0 = 4.2, omega0 = 1, alpha0 = -3.5, 
    xi1 = 5, omega1 = 1.5, alpha1 = 3.5
  )
  
  # create index given the iteration number, e.g., 0001, 0002, ..., 1000 
  setnum <- sprintf("%04d.rds", i)
  
  # save datasets
  saveRDS(df1_sn,
          file = paste0("./data/simdata-setting1/setting1-sn-", setnum))
  saveRDS(df2_sn,
          file = paste0("./data/simdata-setting2/setting2-sn-", setnum))
  saveRDS(df3_sn,
          file = paste0("./data/simdata-setting3/setting3-sn-", setnum))
}
