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
  filepath <- system.file(fname, package = "lgpr2")
  code <- read_file_lines(filepath)
  paste0("  // --------------- ", name, " ---------------\n", code, "\n")
}

# Read lines from file
read_file_lines <- function(file) {
  a <- readLines(file)
  paste(a, collapse = "\n")
}

# Likelihood in model block
stancode_lon_likelihood <- function(stanname_y, dataname) {
  code <- "  // Gaussian likelihood\n"
  paste0(
    code, "  real log_lik_lon = normal_lpdf(", stanname_y, " | f_sum_",
    dataname, ", sigma);"
  )
}

# Log likelihood in generated quantities
stancode_lon_gq <- function(mod, stanname_y, sc_terms_gq) {
  y_var <- mod$y_var
  def_ll <- paste0("  vector[n_LON] log_lik;")
  sylp <- paste0(y_var, "_log_pred")
  def_yp <- paste0("  vector[n_LON] ", sylp, ";")
  line_yp <- paste0(
    "    ", sylp, "[i] = normal_rng(", "f_sum[i], sigma);"
  )
  loop <- paste0("  for(i in 1:n_LON) {\n", line_yp, "\n  }")
  paste(sc_terms_gq, "\n  // Other generated quantities",
    def_ll, def_yp, loop,
    sep = "\n"
  )
}

# Data block
stancode_lon_data <- function(stanname_y, dataname) {
  line0 <- "  // Observation model"
  y_decl <- paste0("  vector<lower=0>[n_", dataname, "] ")
  line1 <- paste0(y_decl, stanname_y, "; // TS model observations")
  paste(line0, line1, sep = "\n")
}

# Transformed data block
stancode_lon_tdata <- function(stanname_y, dataname) {
  c1 <- paste0("  real y_loc = mean(", stanname_y, ");")
  c2 <- paste0("  real y_scale = sd(", stanname_y, ");")
  paste(c1, c2, sep = "\n")
}

# Stan code always in the model block
stancode_loglik <- function(model_suffixes) {
  ll <- paste0("log_lik_", model_suffixes, collapse = " + ")
  code <- "\n  // Likelihood\n"
  paste0(code, "  if(prior_only == 0){\n    target += ", ll, ";\n  }")
}
