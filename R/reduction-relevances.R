compute_noise_amount <- function(y, h) {
  resid <- y - h
  signal_var <- stats::var(h)
  error_var <- sum(resid^2) / (length(h) - 1)
  p_signal <- signal_var / (signal_var + error_var)
  1 - p_signal
}


# Projection
project_draws <- function(model, dat, h_df, formula) {
  checkmate::assert_class(model, "LonModel")
  checkmate::assert_class(dat, "data.frame")
  checkmate::assert_class(h_df, "data.frame")
  checkmate::assert_class(formula, "formula")
  udidx <- unique(h_df$.draw_idx)
  S <- length(udidx)
  pb <- progress::progress_bar$new(total = S)
  h_df[[model$y_var]] <- h_df$value # ref model pred as data
  projs <- list()
  for (s in 1:S) {
    pb$tick()
    h_s <- h_df %>% dplyr::filter(.draw_idx == udidx[s])
    projs[[s]] <- project_draw(formula, h_s)
  }
  projs
}

# Project single draw
project_draw <- function(formula, df) {
  gam_fit <- project_gam.mgcv(formula, df)
  h_proj <- as.numeric(predict(gam_fit))
  h_proj
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
