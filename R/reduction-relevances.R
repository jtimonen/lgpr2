# Signal-to-noise decomposition
compute_noise_amount <- function(y, h) {
  resid <- y - h
  signal_var <- stats::var(h)
  error_var <- sum(resid^2) / (length(h) - 1)
  p_signal <- signal_var / (signal_var + error_var)
  1 - p_signal
}
