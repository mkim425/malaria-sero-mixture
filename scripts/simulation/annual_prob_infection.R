# Annual probability of infection

# Parameters
# lower: the infection rate at the end of the time period
# upper: the infection rate at the beginning of the time period
# beta: a parameter that shifts the curve's transition earlier or later in time
# gamma: a parameter that governs the steepness of the transition
# t0: sigmoid midpoint
annual_prob_infection <- function(lower, upper, beta, gamma, t0) {
  s <- function(t) lower + (upper - lower) / (1 + beta * exp(-gamma * (-t + t0)))
  return(s)
}
