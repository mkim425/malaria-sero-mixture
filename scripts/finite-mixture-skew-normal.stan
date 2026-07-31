/*
  Finite mixture model for Malaria serology data
  With changing probability of infection over time, s(t)
  Use splines for s(t) : logit(s(t)) = c * B
  where c is the weight for the basis functions in B
  Note: RW2 prior on c_i, Numerically stable skew-normal log-PDF, NCRE on RW2
  Individual-level random effect u_i ~ N(0,1) marginalized analytically:
    ever-infected component becomes SN(xi2, omega*, alpha*) where
    omega* = sqrt(omega2^2 + eta^2)
    alpha* = omega2 * alpha2 / sqrt(omega2^2 + eta^2 * (1 + alpha2^2))
  Half-Cauchy(0,0.5) prior on eta
  mean of each skew normal distribution:
  mu0 = mean of SN(location[1], scale[1], shape[1])
  mu1 = mean of SN(location[2], omega*, shape*)
*/

// function to improve numerical stability when alpha*z << 0
functions {
  real safe_skew_normal_lpdf(real y, real xi, real omega, real alpha) {
    real z = (y - xi) / omega;  // standardized variable
    real v_safe = log_sum_exp(alpha * z, -30.0); // smooth lower bound on alpha*z
    return log(2.0) - log(omega) + std_normal_lpdf(z) + normal_lcdf(v_safe | 0, 1);
  }
}

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  vector[N]  y;            // log transformed value of observations
  array[N] int a;          // ages
  int num_basis;           // number of basis functions
  int num_times;           // number of time points
  matrix[num_times, num_basis] B;  // basis matrix
  }


// xi = location (real)
// omega = scale (positive, real)
// alpha = shape (real)

parameters {
  // --- finite mixture (K=1:never-infected/ K=2: ever-infected)
  ordered[K] location;    // locations of mixture components
  vector<lower=0.05>[K] scale;    // scales of mixture components
  vector[K] shape;    // shapes of mixture components
  // --- scale of individual-level random effects
  real<lower=0> eta;    // half-Cauchy(0,1)
  // --- P-spline: non-centred RW2 parameterisation
  // initial conditions of the random walk
  real c1;
  real c2;
  vector[num_basis - 2] z_k;
  real<lower=0> tau;                 // RW2 scale
}


// add transformed_parameters
transformed parameters{
  // --- spline
  vector[num_basis] c;            // spline coefficients on logit scale
  vector<lower=0, upper=1>[num_times] s_hat; // probability of infection
  // --- infection probabilities
  vector[N] logtheta; //log of theta, log P(Z=0): log probability of no infection by age a
  // --- marginalized ever-infected parameters (u_i integrated out analytically)
  real omega_star;   // marginal scale: sqrt(omega2^2 + eta^2)
  real alpha_star;   // marginal shape: omega2 * alpha2 / sqrt(omega2^2 + eta^2*(1+alpha2^2))
  // --- likelihood bookkeeping
  vector[N] log_lik; // log likelihood

  // reconstruct c from initial conditions and scaled innovations.
  c[1] = c1;
  c[2] = c2;
  for (k in 3:num_basis)
    c[k] = 2 * c[k-1] - c[k-2] + tau * z_k[k - 2];
  // --- s(t)
  s_hat = inv_logit(B * c);

  for (n in 1:N) {
    // --- cumulative probability of never being infected
    int t_start;
    t_start = max(1, num_times - a[n]); // prevent invalid indexing
    logtheta[n] = sum(log1m(s_hat[t_start:num_times]));
   }

  // --- marginalized SN parameters for ever-infected component
  // Integrating u_i ~ N(0,1) out of SN(xi2 + eta*u_i, omega2, alpha2) gives
  // SN(xi2, omega*, alpha*) with:
  omega_star = sqrt(square(scale[2]) + square(eta));
  alpha_star = scale[2] * shape[2] / sqrt(square(scale[2]) + square(eta) * (1 + square(shape[2])));

  // Compute the mixture log-likelihood for each observation.
  for (n in 1:N) {
    vector[K] lps;  // log probabilities for each component
    real logtheta_safe;  // numerically guarded version of logtheta[n]

    logtheta_safe = fmin(logtheta[n], -1e-10);

    // --- log pdf
    // safe_skew_normal_lpdf guards log Phi(alpha*z) via log_sum_exp(alpha*z, -30),
    // preventing NaN gradients when an observation falls in the tail of a component.
    lps[1] = logtheta_safe + safe_skew_normal_lpdf(y[n] | location[1], scale[1], shape[1]);
    lps[2] = log1m_exp(logtheta_safe) + safe_skew_normal_lpdf(y[n] | location[2], omega_star, alpha_star);
    // log likelihood
    log_lik[n] = log_sum_exp(lps); //store target for use in LOO IC
    }
}


model {
  // --- priors
  // mixture priors
  location[1] ~ normal(0, 10);
  location[2] ~ normal(0, 10);
  scale ~ student_t(3, 0, 10); // df=3, scale=10, truncated to [0.05, inf)
  shape ~ normal(0, 10);   // center 0, sd = 10 
  // half-Cauchy(0,1) prior on eta
  eta ~ cauchy(0, 2); // half-Cauchy due to lower=0 constraint
  // spline priors
  c1 ~ normal(-3, 0.5) ;   // mean=-3 is between logit(0.03) and logit(0.1) i.e., max s(t) on logit scale
  c2 ~ normal(-3, 0.5);
  z_k ~ std_normal();
  tau ~ normal(0, 0.5);
  // --- likelihood ---
  target += sum(log_lik);
}

// add generated quantities
generated quantities {
  vector[N] y_tilde; // predicted data
  array[N] int<lower=0,upper=1> z; // latent infection status, either 0 or 1
  vector[N] log_prob_noinfection;           // log P(Z_n = 0 | y_n)
  vector[N] prob_noinfection;               // P(Z_n = 0 | y_n)
  real mu0;  // mean of never-infected component SN(location[1], scale[1], shape[1])
  real mu1;  // mean of ever-infected component SN(location[2], omega_star, alpha_star)

  // mean of SN(xi, omega, alpha) = xi + omega * delta * sqrt(2/pi), delta = alpha/sqrt(1+alpha^2)
  mu0 = location[1] + scale[1]  * (shape[1]   / sqrt(1 + square(shape[1])))   * sqrt(2.0 / pi());
  mu1 = location[2] + omega_star * (alpha_star / sqrt(1 + square(alpha_star))) * sqrt(2.0 / pi());

  for (n in 1:N) {
    real logtheta_safe;  // numerically guarded logtheta
    vector[K] lps;
    real p_infection;    // P(Z_n = 1): marginal probability of ever being infected

    logtheta_safe = fmin(logtheta[n], -1e-10);

    // --- posterior predictive draws from marginal SN (u integrated out)
    p_infection = 1 - exp(logtheta_safe);
    z[n] = bernoulli_rng(p_infection); // status of infection
    if (z[n] == 0) {
      y_tilde[n] = skew_normal_rng(location[1], scale[1], shape[1]);
    } else {
      y_tilde[n] = skew_normal_rng(location[2], omega_star, alpha_star);
    }

    // --- posterior probability of never being infected
    lps[1] = logtheta_safe            + safe_skew_normal_lpdf(y[n] | location[1], scale[1], shape[1]);
    lps[2] = log1m_exp(logtheta_safe) + safe_skew_normal_lpdf(y[n] | location[2], omega_star, alpha_star);
    log_prob_noinfection[n] = lps[1] - log_sum_exp(lps); // log P(Z=0 | y_n)
    prob_noinfection[n] = exp(log_prob_noinfection[n]);  // P(Z=0 | y_n)
  }
}


