#' Observation model class (R6 class)
#'
#' @export
Likelihood <- R6::R6Class("Likelihood",
  private = list(
    c_hat = 0
  ),
  public = list(

    #' @description
    #' Get c_hat value
    #'
    #' @param value Value for c_hat.
    get_c_hat = function() {
      private$c_hat
    },

    #' @description
    #' Set c_hat based on y
    #'
    #' @param y Data vector.
    set_c_hat = function(y) {
      value <- mean(y)
      private$c_hat <- value
    }
  )
)

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
    #' The model description as a string
    #' @param y_name Name of y variable.
    string = function(y_name = "y") {
      paste0(
        "  ", y_name, " ~ N(h, sigma^2), where\n",
        "  h = f + c_hat"
      )
    },


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
      c0 <- "  // Observation model"
      y_decl <- paste0("  vector[n_", dataname, "] ")
      c1 <- paste0(y_decl, stanname_y, "; // Longitudinal observations")
      c2 <- "  real c_hat; // fixed offset"
      paste(c0, c1, c2, sep = "\n")
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
      c0 <- "  // Gaussian likelihood"
      c1 <- paste0(
        "  vector[n_LON] h_LON = f_sum_LON + c_hat;"
      )
      c2 <- paste0(
        "  real log_lik_lon = normal_lpdf(", stanname_y, " | h_LON, sigma);"
      )
      paste(c0, c1, c2, sep = "\n")
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
      line_h <- paste0("  vector[n_LON] h = h_LON;")
      line_yp <- paste0(
        "    ", sylp, "[i] = normal_rng(", "h[i], sigma);\n"
      )
      line_ll <- paste0(
        "    log_lik[i] = normal_lpdf(",
        stanname_y, "[i] | h_LON[i], sigma);\n"
      )
      loop <- paste0("  for(i in 1:n_LON) {\n", line_yp, line_ll, "\n  }")
      paste("  // Other generated quantities",
        line_h,
        def_ll, def_yp, loop,
        sep = "\n"
      )
    }
  )
)
