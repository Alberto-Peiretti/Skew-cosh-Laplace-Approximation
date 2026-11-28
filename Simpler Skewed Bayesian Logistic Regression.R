rm(list = ls())
#work_dir <- #Set your working directory here
work_dir <- "C:\\Users\\B00828442\\Documents\\Università\\PhD-ESSEC\\Cut-posteriors and skewed\\Simulazione logistica\\Github"
setwd(work_dir)

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

#   # 2) Find the posterior mode and construct the Laplace approximation
#   data_full <- list(
#     N = n,
#     D = d,
#     x = dataset$X,
#     y = dataset$y,
#     sigma = 2
#   )

#   model <- stan_model("Post-logistic.stan")
#   lap <- rstan::optimizing(model, data = data_full, hessian = TRUE, seed = 1)

#   if (is.null(lap) || !isTRUE(lap$return_code == 0)) {
#     stop("Stan optimization failed for ", suffix, call. = FALSE)
#   }

#   m_lap <- as.numeric(lap$par)
#   precision_lap <- -lap$hessian
#   cov_lap <- chol2inv(chol(precision_lap))

#   # 3) Evaluate the logistic log-posterior and its score at the mode
#   prior_var <- data_full$sigma^2

#   log_post <- function(params) {
#     eta <- drop(X %*% params)
#     sum(plogis((2 * y - 1) * eta, log.p = TRUE)) -
#       sum(params^2) / (2 * prior_var)
#   }

#   eta_mode <- drop(X %*% m_lap)
#   residual <- ifelse(y == 1, plogis(-eta_mode), -plogis(eta_mode))
#   score_lap <- drop(crossprod(X, residual)) - m_lap / prior_var

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
#     score = score_lap,
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

#     # 6) Calculate and save the ESS results
#   res <- multiple_ess_compare(
#     nsim = 100,
#     m = m_lap,
#     Cov = cov_lap,
#     log_post = log_post,
#     logP_theta = logP_theta,
#     N = 1000,
#     radius_sq = mahal_radius_sq_cut,
#     model = sampler_model,
#     debug = TRUE,
#     progress_file = tempfile(
#       pattern = paste0("ESS_progress_", suffix, "_"),
#       tmpdir = getwd(),
#       fileext = ".csv"
#     )
#   )
#   save(res, file = paste0("ESS_lap_cosh_raw_data_nb_c0_", suffix))
#   cat(sprintf("Finished configuration: d=%d, n=%d, seed=%d\n", d, n, seed))
#   flush(stdout())
#   suffix
# }

# n_workers <- 16L  # Set the number of parallel workers based on your system's capabilities
# cluster <- parallel::makePSOCKcluster(min(length(jobs), n_workers))
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
check_ess <- function(res, file) {
  for (method in c("ess_abs_lap", "ess_abs_skewlap", "ess_abs_cosh")) {
    x <- res[[method]]

    if (!is.numeric(x) || length(x) != 100L) {
      stop(sprintf("%s, %s: expected 100 numeric ESS values.",
                   basename(file), method), call. = FALSE)
    }

    bad <- which(!is.finite(x) | x < 1 - 1e-8 | x > 1000 + 1e-8)
    if (length(bad) > 0L) {
      stop(sprintf("%s, %s: invalid ESS at repetition(s) %s.",
                   basename(file), method, paste(bad, collapse = ", ")),
           call. = FALSE)
    }
  }
}

read_ess_data <- function(file, seed) {
  
  env <- new.env()
  objects_loaded <- load(file, envir = env)  
  res_i <- env$res
  check_ess(res_i, file)
  
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
    ) 
  
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

#This block keeps only the files from this code run, avoiding a longer table with runs from other codes in the same directory.
file_index <- file_index %>%
  semi_join(
    settings,
    by = c("n_data" = "n", "dimension" = "d", "seed" = "seed")
  )


file_index <- file_index %>%
  group_by(n_data, dimension, seed) %>%
  slice(1) %>%
  ungroup()

ess_statistics <- function(x) {
  
  x <- as.numeric(x)  
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

library(dplyr)
library(knitr)
library(kableExtra)

df_tex <- ess_gt %>% as.data.frame() %>%
  arrange(Dimension, Seed) %>%
  mutate(
    Dimension = as.character(Dimension),
    Dimension = ifelse(
      grepl("^d\\s*=", Dimension),
      Dimension,
      paste0("d = ", Dimension)
    ),
    n = format(
      as.numeric(gsub(",", "", as.character(n))),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    Gaussian_Min = sprintf("%.1f", as.numeric(Gaussian_Min)),
    SkewLaplace_Min = sprintf("%.1f", as.numeric(SkewLaplace_Min)),
    SkewCosh_Min = sprintf("%.1f", as.numeric(SkewCosh_Min))
  )

tabella_tex <- df_tex %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = "",
    align = c("l", "r", "c", "c", "c", "c"),
    col.names = c(
      "Dimension", "n", "Seed",
      "Minimum", "Minimum", "Minimum"
    )
  ) %>%
  add_header_above(
    c(
      " " = 3,
      "Gaussian Laplace" = 1,
      "Skew-symmetric Laplace" = 1,
      "Skew-cosh Laplace" = 1
    )
  ) %>%
  add_header_above(
    c("Minima ESS of importance weights" = 6)
  ) %>%
  kable_styling(
    font_size = 10,
    full_width = FALSE
  ) %>%
  collapse_rows(
    columns = 1,
    valign = "top",
    latex_hline = "major"
  )

writeLines(
  as.character(tabella_tex),
  "Table_minima_ESS.tex"
)