# Autoformat a 'Stan' code string
autoformat_stancode <- function(code) {
  tryCatch(
    {
      file <- cmdstanr::write_stan_file(code)
      model <- cmdstanr::cmdstan_model(file, compile = FALSE)
      res <- processx::run(
        file.path(cmdstanr::cmdstan_path(), "bin", "stanc"),
        args = c(model$stan_file(), "--auto-format")
      )
      return(res$stdout)
    },
    error = function(e) {
      stop(e$stderr)
    }
  )
}

# Read stan code from inst/functions dir
read_stan_function <- function(name) {
  checkmate::assert_character(name, min.chars = 1)
  fname <- file.path("functions", paste0(name, ".stan"))
  filepath <- system.file(fname, package = "sfgp")
  code <- read_file_lines(filepath)
  paste0("  // --------------- ", name, " ---------------\n", code, "\n")
}

# Read lines from file
read_file_lines <- function(file) {
  a <- readLines(file)
  paste(a, collapse = "\n")
}

# Bounds
stancode_bounds <- function(lower = NULL, upper = NULL) {
  if (is.null(lower) && is.null(upper)) {
    return("")
  }
  if (is.null(lower) && !is.null(upper)) {
    return(paste0("<upper=", upper, ">"))
  }
  if (!is.null(lower) && is.null(upper)) {
    return(paste0("<lower=", lower, ">"))
  }
  if (!is.null(lower) && !is.null(upper)) {
    return(paste0("<lower=", lower, ", upper=", upper, ">"))
  }
}

# Likelihood in model block
stancode_ts_likelihood <- function(stanname_y, dataname) {
  code <- "  // TS model likelihood\n"
  paste0(
    code, "  real log_lik_lon = normal_lpdf(", stanname_y, "_log | f_sum_",
    dataname, ", sigma);"
  )
}

# Log likelihood in generated quantities
stancode_ts_gq <- function(mod, stanname_y, sc_terms_gq) {
  cap_val <- mod$log_y_cap
  y_var <- mod$y_var
  def_ll <- paste0("  vector[n_LON] log_lik;")
  sylp <- paste0(y_var, "_log_pred")
  def_yp <- paste0("  vector[n_LON] ", sylp, ";")
  def_fcap <- paste0(
    "  vector[n_LON] f_sum_capped = fmin(f_sum, ", cap_val, ");"
  )
  def_ycap <- paste0("  vector[n_LON] ", y_var, "_log_pred_capped;")
  line_yp <- paste0(
    "    ", sylp, "[i] = normal_rng(", "f_sum[i], sigma);"
  )
  loop <- paste0("  for(i in 1:n_LON) {\n", line_yp, "\n  }")
  line_ycap <- paste0(
    "  ", y_var, "_log_pred_capped = fmin(", sylp, ", ",
    cap_val, ");"
  )
  paste(sc_terms_gq, "\n  // Other generated quantities",
    def_ll, def_yp,
    def_fcap, def_ycap, loop, line_ycap,
    sep = "\n"
  )
}

# Data block
stancode_ts_data <- function(stanname_y, dataname) {
  line0 <- "  // Observation model"
  line1 <- "  real<lower=0> delta; // Offset for log transform of y"
  y_decl <- paste0("  vector<lower=0>[n_", dataname, "] ")
  line2 <- paste0(y_decl, stanname_y, "; // TS model observations")
  paste(line0, line1, line2, sep = "\n")
}

# Transformed data block
stancode_ts_tdata <- function(stanname_y, dataname) {
  c1 <- paste0(
    "  vector[n_", dataname, "] ", stanname_y, "_log = log(", stanname_y,
    " + delta); // log-scale TS observations\n"
  )
  c2 <- paste0("  real y_loc = mean(", stanname_y, ");")
  c3 <- paste0("  real y_scale = sd(", stanname_y, ");")
  paste(c1, c2, c3, sep = "\n")
}

# Stan code always in the model block
stancode_loglik <- function(model_suffixes) {
  ll <- paste0("log_lik_", model_suffixes, collapse = " + ")
  code <- "\n  // Likelihood\n"
  paste0(code, "  if(prior_only == 0){\n    target += ", ll, ";\n  }")
}
