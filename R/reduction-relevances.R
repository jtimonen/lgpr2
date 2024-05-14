# Signal-to-noise decomposition
compute_noise_amount <- function(y, h) {
  resid <- y - h
  signal_var <- stats::var(h)
  error_var <- sum(resid^2) / (length(h) - 1)
  p_signal <- signal_var / (signal_var + error_var)
  1 - p_signal
}


# Projection
project_draws <- function(fit, dat, h_df, formula) {
  checkmate::assert_class(fit, "LonModelFit")
  checkmate::assert_class(dat, "data.frame")
  checkmate::assert_class(h_df, "data.frame")
  checkmate::assert_class(formula, "formula")

  # Extract things
  model <- fit$get_model()
  udidx <- unique(h_df$.draw_idx)
  sigma_ref <- posterior::as_draws_array(fit$draws("sigma"))
  y_data <- dat[[model$stanname_y()]] # actual data

  # Loop over draws
  S <- length(udidx)
  pb <- progress::progress_bar$new(total = S)
  projs <- list()
  for (s in 1:S) {
    dat_mgcv <- dat
    pb$tick()
    h_s <- h_df %>% dplyr::filter(.draw_idx == udidx[s])
    mu_ref <- h_s$value
    dat_mgcv[[model$stanname_y()]] <- mu_ref
    projs[[s]] <- project_draw(formula, y_data, mu_ref, dat_mgcv, sigma_ref[s])
  }
  list(
    gam_fits = sapply(projs, function(x) x$gam_fit),
    mu_ref = sapply(projs, function(x) x$mu_ref),
    mu_proj = sapply(projs, function(x) x$mu_proj),
    loglik_ref = sapply(projs, function(x) x$loglik_ref),
    loglik_proj = sapply(projs, function(x) x$loglik_proj),
    kl_div = sapply(projs, function(x) x$kl_div)
  )
}

# Project single draw
project_draw <- function(formula, y_dat, mu_ref, df, sigma_ref) {
  gam_fit <- project_gam.mgcv(formula, df)
  mu_proj <- as.numeric(predict(gam_fit))
  sigma_proj <- project_sigma(sigma_ref, mu_ref, mu_proj)
  kl_div <- kl_divergence_gaussians(mu_ref, mu_proj, sigma_ref, sigma_proj)
  loglik_ref <- loglik_gaussian(y_dat, mu_ref, sigma_ref)
  loglik_proj <- loglik_gaussian(y_dat, mu_proj, sigma_proj)

  # Return
  list(
    gam_fit = gam_fit,
    kl_div = kl_div,
    mu_proj = mu_proj,
    mu_ref = mu_ref,
    loglik_ref = loglik_ref,
    loglik_proj = loglik_proj
  )
}


# Project to a GAM using mgcv
project_gam.mgcv <- function(form, dat, ...) {
  min_sp <- getOption("min.sp", default = NULL)
  if (!is.null(min_sp)) {
    D <- length(attr(terms(form), "term.labels"))
    min_sp <- rep(min_sp, D)
  }
  mgcv::gam(form, data = dat, min.sp = min_sp, ...)
}


# Euclidean norm squared
squared_euc_norm <- function(x) {
  v <- as.vector(x)
  sum(v * v)
}

# Project sigma
project_sigma <- function(sigma, mu, mu_proj) {
  v <- mu_proj - mu
  n <- length(v)
  sigma_proj2 <- sigma^2 + 1 / n * squared_euc_norm(v)
  sqrt(sigma_proj2)
}

# KL divergence of multivariate normals with S = s^2 I
kl_divergence_gaussians <- function(mu, mu_proj, sigma, sigma_proj) {
  term1 <- log(sigma_proj / sigma) + (sigma^2) / (2 * sigma_proj^2) - 0.5
  term2 <- squared_euc_norm(mu - mu_proj) / (2 * sigma_proj^2)
  length(mu) * term1 + term2
}

# log likelihood at each point
loglik_gaussian <- function(y, mu, sigma) {
  stats::dnorm(y,
    mean = mu,
    sd = sigma,
    log = TRUE
  )
}
