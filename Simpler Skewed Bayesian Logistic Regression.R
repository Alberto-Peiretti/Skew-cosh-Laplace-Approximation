rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(gt)
library(reshape)

rstan::rstan_options(auto_write = TRUE)
invisible(rstan::stan_model("Post-logistic.stan"))
t0 <- proc.time()
configurations <- data.frame(
  d = c(20L, 30L, 40L, 50L),
  n = c(2100L, 4400L, 9500L, 12500L)
)
seeds <- c(1L, 3L, 4L, 7L)
settings <- expand.grid(
  config = seq_len(nrow(configurations)),
  seed = seeds,
  KEEP.OUT.ATTRS = FALSE
)
settings$d <- configurations$d[settings$config]
settings$n <- configurations$n[settings$config]
jobs <- lapply(seq_len(nrow(settings)), function(i) {
  as.list(settings[i, c("d", "n", "seed")])
})

# run_simulation <- function(job) {
#   d <- job$d
#   n <- job$n
#   seed <- job$seed
#   suffix <- sprintf("n%d_d%d_%d", n, d, seed)
#   cat(sprintf("Starting configuration: d=%d, n=%d, seed=%d\n", d, n, seed))
#   flush(stdout())

#   # 1) Simulate logistic regression data with normally distributed coefficients
#   set.seed(seed)
#   X <- matrix(runif(n * d, min = 0, max = 1), nrow = n, ncol = d)
#   beta <- rnorm(d, mean = 0, sd = 1)
#   alpha <- 0
#   eta <- as.numeric(alpha + X %*% beta)
#   p_vec <- plogis(eta)
#   y <- rbinom(n, size = 1, prob = p_vec)
#   dataset <- list(d = d, n = n, X = X, y = y)

#   # 2) Calculate the posterior
#   data_full <- list(
#     N = n,
#     D = d,
#     x = dataset$X,
#     y = dataset$y,
#     sigma = 2
#   )
#   post <- stan(
#     file = "Post-logistic.stan",
#     data = data_full,
#     warmup = 1000,
#     iter = 10^4,
#     chains = 4,
#     seed = 1
#   )
#   save(post, file = paste0("Posterior_Logistic_HMC_", suffix))

#   # 3) Construct the Laplace approximation
#   model <- stan_model("Post-logistic.stan")
#   lap <- rstan::optimizing(model, data = data_full, hessian = TRUE, seed = 1)
#   m_lap <- lap$theta_tilde[1:d]
#   cov_lap <- -solve(lap$hessian)

#   # 4) Define the log-posterior function
#   fit <- rstan::vb(
#     model,
#     data = data_full,
#     iter = 10000,
#     algorithm = "fullrank",
#     tol_rel_obj = 1e-8,
#     seed = 1,
#     output_samples = 10^4,
#     adapt_iter = 300
#   )
#   log_post <- function(params) {
#     rstan::log_prob(fit, params, adjust_transform = TRUE)
#   }

#   # 5) Construct the skew-cosh envelope used by the sampler
#   logcosh_stable <- function(x) {
#     ax <- abs(x)
#     ax + log1p(exp(-2 * ax)) - log(2)
#   }
#   mahal_radius_sq_cut <- qchisq(1 - 1e-3, df = d)
#   inside_ball <- function(theta, center = m_lap, Sigma = cov_lap,
#                           radius_sq = mahal_radius_sq_cut) {
#     mahalanobis(theta, center = center, cov = Sigma) < radius_sq
#   }
#   logP_theta <- function(theta, center = m_lap, Sigma = cov_lap,
#                          radius_sq = mahal_radius_sq_cut) {
#     if (!inside_ball(theta, center = center, Sigma = Sigma,
#                      radius_sq = radius_sq)) {
#       return(0)
#     }
#     theta_ref <- 2 * center - theta
#     delta <- 0.5 * (log_post(theta) - log_post(theta_ref))
#     logcosh_stable(delta)
#   }

#   sampler_model <- make_skew_cosh(
#     X = X,
#     y = y,
#     center = m_lap,
#     Sigma = cov_lap,
#     radius_sq = mahal_radius_sq_cut,
#     log_post = log_post,
#     score = as.numeric(rstan::grad_log_prob(
#       fit, m_lap, adjust_transform = TRUE
#     )),
#     quadratic_prior = TRUE,
#     partition = "radial",
#     search = "cosh"
#   )
#   sampler_model <- make_skew_cosh_directional(
#     model = sampler_model,
#     n_radial_splits = 32L,
#     max_steps = 16L,
#     angular_bins = 128L,
#     integrated_bound = TRUE,
#     quadratic_remainder = TRUE,
#     min_log_bound_improvement = 0.001,
#     stagnation_patience = 1L,
#     verbose = TRUE
#   )
#   sampler_model <- refine_skew_cosh_directional_local(
#     model = sampler_model,
#     max_cells = Inf,
#     n_radial_steps = 32L,
#     angular_bins = 128L,
#     local_bins = 128L,
#     min_log_bound_improvement = 0.001,
#     stagnation_patience = 1L,
#     verbose = TRUE
#   )
#   sampler_model <- refine_skew_cosh_directional_angular(
#     model = sampler_model,
#     max_cells = Inf,
#     n_radial_steps = 32L,
#     angular_bins = 128L,
#     local_bins = 128L,
#     min_log_bound_improvement = 0.001,
#     stagnation_patience = 1L,
#     verbose = TRUE
#   )

#   # 6) Calculate and save the ESS results
#   res <- multiple_ess_compare(
#     nsim = 100,
#     m = m_lap,
#     Cov = cov_lap,
#     log_post = log_post,
#     logP_theta = logP_theta,
#     post = post,
#     N = 1000,
#     radius_sq = mahal_radius_sq_cut,
#     model = sampler_model,
#     debug = TRUE
#   )
#   save(res, file = paste0("ESS_lap_cosh_raw_data_nb_c0_", suffix))
#   cat(sprintf("Finished configuration: d=%d, n=%d, seed=%d\n", d, n, seed))
#   flush(stdout())
#   suffix
# }

# cluster <- parallel::makePSOCKcluster(4L)
# results <- tryCatch({
#   parallel::clusterCall(cluster, setwd, work_dir)
#   parallel::clusterEvalQ(cluster, {
#     source("Functions.R")
#     source("skew_cosh.R")
#     library(StanHeaders)
#     library(rstan)
#     options(mc.cores = 4L)
#     rstan_options(auto_write = TRUE)
#     NULL
#   })
#   parallel::parLapplyLB(cluster, jobs, run_simulation, chunk.size = 1L)
# }, finally = parallel::stopCluster(cluster))
# print((proc.time() - t0)["elapsed"])

# Plot ESS importance weights different datasets (Figure 2-3-4-5)

read_ess_data <- function(file, seed) {
  
  env <- new.env()
  objects_loaded <- load(file, envir = env)  
  res_i <- env$res
  
  data.frame(
    Gaussian_Laplace = res_i$ess_abs_lap,
    Skew_Laplace     = res_i$ess_abs_skewlap,
    Skew_Cosh        = res_i$ess_abs_cosh
  ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Method",
      values_to = "ESS"
    ) %>%
    mutate(
      Seed = seed
    )
}


make_ess_plot <- function(d, n) {
  
  files <- setNames(
    paste0(
      "ESS_lap_cosh_raw_data_nb_c0_n",
      n,
      "_d",
      d,
      "_",
      seeds
    ),
    seeds
  )
  
  df_hist_all <- bind_rows(
    Map(
      f = read_ess_data,
      file = unname(files),
      seed = names(files)
    )
  ) %>%
    mutate(
      Method = factor(
        Method,
        levels = c(
          "Gaussian_Laplace",
          "Skew_Laplace",
          "Skew_Cosh"
        ),
        labels = c(
          "Gaussian Laplace",
          "Skew-symmetric Laplace",
          "Skew-cosh Laplace"
        )
      ),
      Seed = factor(
        Seed,
        levels = names(files),
        labels = paste0("Seed = ", names(files))
      )
    ) %>%
    filter(is.finite(ESS))
  
  p_ess <- ggplot(
    df_hist_all,
    aes(
      x = ESS,
      fill = Method,
      color = Method
    )
  ) +
    geom_histogram(
      bins = 50,
      alpha = 0.35,
      position = "identity",
      drop = "all"
    ) +
    facet_wrap(
      ~ Seed,
      nrow = 2,
      ncol = 2,
      scales = "free"
    ) +
    scale_fill_manual(
      values = c(
        "Gaussian Laplace" = "#E41A1C",
        "Skew-symmetric Laplace" = "#377EB8",
        "Skew-cosh Laplace" = "#4DAF4A"
      )
    ) +
    scale_color_manual(
      values = c(
        "Gaussian Laplace" = "#E41A1C",
        "Skew-symmetric Laplace" = "#377EB8",
        "Skew-cosh Laplace" = "#4DAF4A"
      )
    ) +
    theme_minimal(base_size = 22) +
    theme(
    legend.position = "bottom",
    legend.title = element_text(size = 22, face = "plain"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "plain"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "plain"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.spacing = grid::unit(1, "lines")
    ) +
    labs(
      x = "ESS",
      y = "Frequency",
      title = paste0(
        "ESS importance weights: d = ",
        d,
        ", n = ",
        n
      ),
      fill = "Approximation",
      color = "Approximation"
    )
  
  print(p_ess)
  
  ggsave(
  filename = paste0(
    "ESS_weights_d",
    d,
    "_n",
    n,
    ".pdf"
  ),
  plot = p_ess,
  width = 16,
  height = 9,
  units = "in"
  )
}

# Generate and save plots for all configurations
for (i in seq_len(nrow(configurations))) {
  make_ess_plot(
    d = configurations$d[i],
    n = configurations$n[i]
  )
}

#ESS minima importance weights table (Table 2)
files_found <- list.files(
  path = work_dir,
  pattern = "^ESS_lap_cosh_raw_data_nb_c0_n",
  full.names = TRUE
)

parsed_names <- str_match(
  basename(files_found),
  "_n([0-9]+)_d([0-9]+)_([0-9]+)(?:\\([0-9]+\\))?(?:\\.[^.]+)?$"
)

file_index <- tibble(
  file = files_found,
  filename = basename(files_found),
  n_data = as.integer(parsed_names[, 2]),
  dimension = as.integer(parsed_names[, 3]),
  seed = as.integer(parsed_names[, 4])
) %>%
  filter(
    !is.na(n_data),
    !is.na(dimension),
    !is.na(seed)
  )

file_index <- file_index %>%
  mutate(
    is_copy = str_detect(
      filename,
      "\\([0-9]+\\)(?:\\.[^.]+)?$"
    )
  ) %>%
  arrange(
    dimension,
    n_data,
    seed,
    is_copy,
    filename
  )

file_index <- file_index %>%
  group_by(n_data, dimension, seed) %>%
  slice(1) %>%
  ungroup()

ess_statistics <- function(x) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(
      c(
        min = NA_real_
      )
    )
  }
  
  c(
    min = min(x)
  )
}

read_ess_summary <- function(
    file,
    n_data,
    dimension,
    seed
) {
  
  env <- new.env()
  loaded_objects <- load(file, envir = env)
  
  res_i <- env$res
  
  required_components <- c(
    "ess_abs_lap",
    "ess_abs_skewlap",
    "ess_abs_cosh"
  )
  
  stat_lap <- ess_statistics(
    res_i$ess_abs_lap
  )
  
  stat_skewlap <- ess_statistics(
    res_i$ess_abs_skewlap
  )
  
  stat_cosh <- ess_statistics(
    res_i$ess_abs_cosh
  )
  
  tibble(
    n = n_data,
    d = dimension,
    seed = seed,
    
    Gaussian_min = stat_lap["min"],
    SkewLaplace_min = stat_skewlap["min"],
    SkewCosh_min = stat_cosh["min"]
  )
}

ess_summary_numeric <- pmap_dfr(
  file_index %>%
    select(
      file,
      n_data,
      dimension,
      seed
    ),
  read_ess_summary
) %>%
  arrange(d, seed, n)

ess_summary_table <- ess_summary_numeric %>%
  transmute(
    Dimension = paste0("d = ", d),
    n = n,
    Seed = seed,
    
    Gaussian_Min = Gaussian_min,
    SkewLaplace_Min = SkewLaplace_min,
    SkewCosh_Min = SkewCosh_min
  )

ess_gt <- ess_summary_table %>%
  gt(
    groupname_col = "Dimension"
  ) %>%
  tab_header(
    title = "Minima ESS of importance weights"
  ) %>%
  tab_spanner(
    label = "Gaussian Laplace",
    columns = c(
      Gaussian_Min
    )
  ) %>%
  tab_spanner(
    label = "Skew-symmetric Laplace",
    columns = c(
      SkewLaplace_Min
    )
  ) %>%
  tab_spanner(
    label = "Skew-cosh Laplace",
    columns = c(
      SkewCosh_Min
    )
  ) %>%
  cols_label(
    n = "n",
    Seed = "Seed",
    Gaussian_Min = "Minimum",
    SkewLaplace_Min = "Minimum",
    SkewCosh_Min = "Minimum"
  ) %>%
  fmt_number(
    columns = c(
      Gaussian_Min,
      SkewLaplace_Min,
      SkewCosh_Min
    ),
    decimals = 1
  ) %>%
  fmt_integer(
    columns = c(n, Seed)
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_stubhead(
    label = "Dimension"
  ) %>%
  tab_options(
    row_group.as_column = TRUE,
    row_group.font.weight = "bold",
    table.font.size = px(13)
  )

show(ess_gt)

