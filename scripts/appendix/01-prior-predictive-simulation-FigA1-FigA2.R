# Prior predictive simulation
#
# Model summary:
#   - Two-component skew-normal mixture: never-infected (K=1) / ever-infected (K=2)
#   - s(t): annual P(infection) modeled with P-spline + non-centred RW2 prior on coefficients
#   - Individual random effect u_i ~ N(0,1) marginalized analytically:
#       ever-infected component -> SN(xi2, omega*, alpha*)
#       omega* = sqrt(omega2^2 + eta^2)
#       alpha* = omega2 * alpha2 / sqrt(omega2^2 + eta^2 * (1 + alpha2^2))
#   - eta ~ half-Cauchy(0, 2)

library(tidyverse)
library(splines)
library(sn)       # rsn() for skew-normal draws
library(gridExtra)

set.seed(345)

# Step 1. Data structure 
n_prior_draws <- 500       # number of prior parameter draws
n_obs         <- 300       # synthetic observations per draw
ages          <- sample(1:95, n_obs, replace = TRUE)  # ages in years

time_points <- 1:100       # calendar years (year 100 = survey year)
num_times   <- length(time_points)
B           <- bs(time_points, df = 6, intercept = TRUE)   # B-spline basis
num_basis   <- ncol(B)

# Step 2. Prior samplers 
# random draw from student_t(df, 0, sigma) truncated at lower=0.05
rtrunc_t <- function(n, df, sigma = 1, lower = 0.05) {
  out <- numeric(n) # pre-allocate output vector
  i   <- 0L         # count of accepted draws
  while (i < n) { 
    x <- sigma * rt(1, df) # draw from student_t(0, sigma) with df degrees of freedom
    if (x >= lower) { i <- i + 1L; out[i] <- x } # accept if above lower bound
  }
  out
}

# random draw from half-Normal: N(0, sd) truncated at 0 
rhalf_normal <- function(n, sd = 1) abs(rnorm(n, 0, sd))

# random draw from half-Cauchy: Cauchy(0, scale) truncated at 0 
rhalf_cauchy <- function(n, scale = 1) abs(rcauchy(n, 0, scale))

# Step 3. One prior predictive draw (n_obs synthetic observation)
draw_prior_predictive <- function(ages, B, num_basis, num_times) {
  n <- length(ages)
  # --- Mixture component parameters ---
  # ordered[K]: location[1] < location[2]
  repeat {
    loc <- c(rnorm(1, 0, 10), rnorm(1, 0, 10))
    if (loc[1] < loc[2]) break
  }
  scale <- rtrunc_t(2, df = 3, sigma = 10, lower = 0.05)   # student_t(3,0,10) | [0.05, inf)
  shape <- rnorm(2, 0, 10)                                   # N(0, 10)

  # --- half-Cauchy(0, 2) on eta (scale of individual random effect) ---
  eta <- rhalf_cauchy(1, scale = 2)

  # --- P-spline: non-centred RW2 reconstruction ---
  c1  <- rnorm(1, -3, 0.5)                     # N(-3, 0.5)
  c2  <- rnorm(1, -3, 0.5)
  z_k <- rnorm(num_basis - 2, 0, 1)            # std_normal innovations
  tau <- rhalf_normal(1, sd = 0.5)             # half-N(0, 0.5)

  c <- numeric(num_basis)    # pre-allocate RW2 coefficients
  c[1] <- c1; c[2] <- c2     # initialize first two coefficients
  for (k in 3:num_basis)     # construct c[k] from non-centred RW2
    c[k] <- 2 * c[k-1] - c[k-2] + tau * z_k[k - 2]

  s_hat <- plogis(as.vector(B %*% c))   # annual P(infection); length = num_times

  # log P(never infected by age a_i): sum log(1 - s(t)) over years alive
  logtheta <- vapply(ages, function(a) {
    t_start <- max(1L, num_times - a)
    sum(log1p(-s_hat[t_start:num_times]))
  }, numeric(1))

  # --- Marginalized ever-infected SN parameters (u_i integrated out) ---
  omega_star <- sqrt(scale[2]^2 + eta^2)
  alpha_star <- scale[2] * shape[2] / sqrt(scale[2]^2 + eta^2 * (1 + shape[2]^2))

  # --- Simulate observations ---
  p_inf <- 1 - exp(logtheta)                  # P(ever infected) by age
  z     <- rbinom(n, 1, p_inf)               # latent infection status

  y <- numeric(n)
  never_idx  <- which(z == 0)
  ever_idx   <- which(z == 1)
  if (length(never_idx) > 0)
    y[never_idx] <- rsn(length(never_idx),
                        xi = loc[1], omega = scale[1],   alpha = shape[1])
  if (length(ever_idx) > 0)
    y[ever_idx]  <- rsn(length(ever_idx),
                        xi = loc[2], omega = omega_star, alpha = alpha_star)

  list(
    y          = y,
    z          = z,
    p_inf      = p_inf,
    s_hat      = s_hat,
    loc        = loc,
    scale      = scale,
    shape      = shape,
    eta        = eta,
    omega_star = omega_star,
    alpha_star = alpha_star,
    tau        = tau,
    c          = c
  )
}

# Step 4. Run simulation 
prior_draws <- replicate(
  n_prior_draws,
  draw_prior_predictive(ages, B, num_basis, num_times),
  simplify = FALSE
)

# Step 5. Collect results
# All simulated y values (prior predictive distribution of the observable)
y_df <- map_dfr(seq_along(prior_draws), function(i) {
  tibble(draw = i, 
         y = prior_draws[[i]]$y,
         infected = factor(prior_draws[[i]]$z,
                           levels = 0:1, labels = c("never", "ever")))
})

# s(t) curves over time
s_df <- map_dfr(seq_along(prior_draws), function(i) {
  tibble(draw = i, time = time_points, s_hat = prior_draws[[i]]$s_hat)
})

# Scalar parameter draws
param_df <- map_dfr(prior_draws, function(d) {
  tibble(
    loc1       = d$loc[1],
    loc2       = d$loc[2],
    scale1     = d$scale[1],
    scale2     = d$scale[2],
    shape1     = d$shape[1],
    shape2     = d$shape[2],
    eta        = d$eta,
    omega_star = d$omega_star,
    alpha_star = d$alpha_star,
    tau        = d$tau
  )
})


# Step  6.
# Prior predictive distribution of y (log antibody titer)
make_y_plot <- function(i) {
  y_df_single_draw <- y_df |> filter(draw == i)
  qs <- quantile(y_df_single_draw$y, c(0.025, 0.975), na.rm = TRUE)
  ggplot(y_df_single_draw, aes(x = y)) +
    geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
    labs(title = paste0("Draw ", i)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 12)
    )
}

p_y <- gridExtra::arrangeGrob(
  grobs = lapply(1:10, make_y_plot),
  ncol = 2,
  top = grid::textGrob("Prior predictive log antibody titer", 
                       gp = grid::gpar(fontsize = 14, fontface = "bold")),
  bottom = grid::textGrob("log(y)", 
                          gp = grid::gpar(fontsize = 12, fontface = "bold")),
  left = grid::textGrob("Density", rot = 90, 
                        gp = grid::gpar(fontsize = 12, fontface = "bold"))
)

# s(t) curves — 95% envelope + median
s_summary <- s_df |>
  group_by(time) |>
  summarise(
    med = median(s_hat),
    mean = mean(s_hat),
    lo  = quantile(s_hat, 0.025),
    hi  = quantile(s_hat, 0.975),
    .groups = "drop"
  )

p_s <- ggplot() +
  geom_line(
    data    = s_df,
    mapping = aes(x = time, y = s_hat, group = draw),
    alpha = 0.18, colour = "steelblue", linewidth = 0.4
  ) +
  geom_ribbon(
    data    = s_summary,
    mapping = aes(x = time, ymin = lo, ymax = hi),
    fill = "steelblue", alpha = 0.25
  ) +
  labs(
    title = "Prior predictive annual probability of infection",
    x = "Year", y = "s(t)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12)
  )


# Step 7. Save plots 
ggsave("plots/FigA1-prior-pred-y-draws-1-10.pdf", p_y, width = 14, height = 18)
ggsave("plots/FigA2-prior-pred-s-curves-full_range.pdf", p_s, width = 7, height = 4.5)




