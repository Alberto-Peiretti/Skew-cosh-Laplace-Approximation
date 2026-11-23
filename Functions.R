# ESS of importance weights for Gaussian Laplace, skew-symmetric Laplace,
# and skew-cosh Laplace. Package required here: mvtnorm.
# Load skew_cosh.R and construct the sampler before calling this function.

w_joint <- function(sim, xi, log_post, ...)
{
  w1 <- log_post(sim,...)
  w2 <- log_post(2*xi -sim,...)
  pr <- 1/(1 + exp( w2 - w1 ))
  pr
}

sim_skew <- function(sim, m, weights)
{
  bern <- rbinom(dim(sim)[1], 1, weights  )
  sim[bern == 0,] <- t( 2 * m - t(sim[bern == 0, ] ) )
  weights[bern == 0] <- 1 - weights[bern == 0]
  return( list( sim = sim, weights = weights, bern = bern  ) )
}

log_phi_lap <- function(x, m, Cov) {
  mvtnorm::dmvnorm(x, mean = m, sigma = Cov, log = TRUE)
}

logiw_laplace <- function(sim, m, Cov, log_post, ...) {
  logpi <- apply(sim, 1, function(z) log_post(z, ...))
  logq  <- log_phi_lap(sim, m, Cov)
  logpi - logq
}

logiw_skew_laplace <- function(sim_skew, w_at_sim_skew, m, Cov, log_post, ...) {
  logpi <- apply(sim_skew, 1, function(z) log_post(z, ...))
  logq  <- log(2) + log_phi_lap(sim_skew, m, Cov) + log(w_at_sim_skew)
  logpi - logq
}

logiw_skew_cosh <- function(sim_cosh, m, Cov, log_post, logP_theta, ...) {
  logpi <- apply(sim_cosh, 1, function(z) log_post(z, ...))
  logP  <- apply(sim_cosh, 1, function(z) logP_theta(z, center = m))
  wj    <- apply(sim_cosh, 1, function(z) w_joint(z, m, log_post, ...))
  
  logq_unnorm <- log(2) + log_phi_lap(sim_cosh, m, Cov) + logP + log(wj)
  logpi - logq_unnorm
}

ess_from_logweights <- function(logw) {
  lw <- as.numeric(logw)
  lw <- lw - max(lw)
  w  <- exp(lw)
  (sum(w)^2) / sum(w^2)
}

multiple_ess_compare <- function(nsim,
                                 m,
                                 Cov,
                                 log_post,
                                 logP_theta,
                                 post = NULL,
                                 N,
                                 radius_sq = model$radius_sq,
                                 debug = FALSE,
                                 debug_iter = integer(0),
                                 model,
                                 max_evals_per_draw = 50000L,
                                 max_seconds_per_draw = 120,
                                 ...) {
  ess_abs_lap <- ess_abs_skewlap <- ess_abs_cosh <- rep(NA_real_, nsim)

  logP_for_weights <- function(theta, ...) {
    logP_theta(theta, center = m, Sigma = Cov, radius_sq = radius_sq)
  }

  for (k in seq_len(nsim)) {
    set.seed(1000 + k)

    # Gaussian Laplace.
    sim_lap <- mvtnorm::rmvnorm(N, mean = m, sigma = Cov)
    logiw_lap <- logiw_laplace(
      sim = sim_lap, m = m, Cov = Cov, log_post = log_post, ...
    )
    ess_abs_lap[k] <- ess_from_logweights(logiw_lap)

    # Skew-symmetric Laplace, using the same Gaussian draws.
    w0 <- apply(sim_lap, 1, function(z) w_joint(z, m, log_post, ...))
    skew_obj <- sim_skew(sim_lap, m, w0)
    logiw_skewlap <- logiw_skew_laplace(
      sim_skew = skew_obj$sim,
      w_at_sim_skew = skew_obj$weights,
      m = m, Cov = Cov, log_post = log_post, ...
    )
    ess_abs_skewlap[k] <- ess_from_logweights(logiw_skewlap)

    # Skew-cosh Laplace.
    cosh_obj <- sample_skew_cosh_rejection(
      N_keep = N, model = model,
      max_evals_per_draw = max_evals_per_draw,
      max_seconds_per_draw = max_seconds_per_draw,
      verbose = FALSE
    )
    logiw_cosh <- logiw_skew_cosh(
      sim_cosh = cosh_obj$sim_skew, m = m, Cov = Cov,
      log_post = log_post, logP_theta = logP_for_weights, ...
    )
    ess_abs_cosh[k] <- ess_from_logweights(logiw_cosh)

    if (debug && (length(debug_iter) == 0L || k %in% debug_iter)) {
      message(sprintf("Ripetizione %d/%d completata.", k, nsim))
    }
  }

  list(
    ess_abs_lap = ess_abs_lap,
    ess_abs_skewlap = ess_abs_skewlap,
    ess_abs_cosh = ess_abs_cosh,
    ess_rel_lap = ess_abs_lap / N,
    ess_rel_skewlap = ess_abs_skewlap / N,
    ess_rel_cosh = ess_abs_cosh / N,
    gain_cosh_vs_lap = 100 * (ess_abs_cosh / ess_abs_lap - 1),
    gain_cosh_vs_skewlap = 100 * (ess_abs_cosh / ess_abs_skewlap - 1),
    gain_skewlap_vs_lap = 100 * (ess_abs_skewlap / ess_abs_lap - 1)
  )
}
