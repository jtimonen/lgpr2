#' Observation model class (R6 class)
#'
#' @export
Likelihood <- R6::R6Class("Likelihood")

#' Observation model class (R6 class)
#'
#' @export
#' @field prior_sigma Prior for noise magnitude parameter.
GaussianLikelihood <- R6::R6Class("GaussianLikelihood",
  inherit = Likelihood,

  # PUBLIC
  public = list(
    prior_sigma = "lognormal(1, 1)",

    #' @description
    #' Create model
    #'
    #' @param prior_sigma Prior for noise magnitude parameter.
    initialize = function(prior_sigma = "lognormal(1, 1)") {
      checkmate::assert_character(prior_sigma)
      self$prior_sigma <- prior_sigma
    },

    #' @description
    #' Stan code for data block
    #'
    #' @param stanname_y Name of y variable in Stan code
    #' @param dataname suffix
    stancode_data = function(stanname_y, dataname) {
      line0 <- "  // Observation model"
      y_decl <- paste0("  vector<lower=0>[n_", dataname, "] ")
      line1 <- paste0(y_decl, stanname_y, "; // Longitudinal observations")
      paste(line0, line1, sep = "\n")
    },

    #' @description
    #' Stan code for parameters block
    stancode_pars = function() {
      paste0("  real<lower=0> sigma; // noise std parameter \n")
    },

    #' @description
    #' Stan code for transformed parameters block
    #' @param stanname_y Name of y variable in Stan code
    #' @param dataname suffix
    stancode_tpars = function(stanname_y, dataname) {
      code <- "  // Gaussian likelihood\n"
      paste0(
        code, "  real log_lik_lon = normal_lpdf(", stanname_y, " | f_sum_",
        dataname, ", sigma);"
      )
    },

    #' @description
    #' Stan code for model block
    #' @param stanname_y Name of y variable in Stan code
    #' @param dataname suffix
    stancode_model = function() {
      c1 <- paste0("  sigma ~ ", self$prior_sigma, ";\n")
      ll <- "log_lik_lon"
      code <- "\n  // Likelihood\n"
      c2 <- paste0(code, "  if(prior_only == 0){\n    target += ", ll, ";\n  }")
      paste(c1, c2, sep = "\n")
    },

    #' @description
    #' Stan code for the generated quantities block
    #' @param stanname_y Name of y variable in Stan code
    #' @param y_name name of y variable
    stancode_gq = function(y_name, stanname_y) {
      y_var <- y_name
      def_ll <- paste0("  vector[n_LON] log_lik;")
      sylp <- paste0(y_var, "_log_pred")
      def_yp <- paste0("  vector[n_LON] ", sylp, ";")
      line_yp <- paste0(
        "    ", sylp, "[i] = normal_rng(", "f_sum[i], sigma);"
      )
      loop <- paste0("  for(i in 1:n_LON) {\n", line_yp, "\n  }")
      paste("  // Other generated quantities",
        def_ll, def_yp, loop,
        sep = "\n"
      )
    }
  )
)
