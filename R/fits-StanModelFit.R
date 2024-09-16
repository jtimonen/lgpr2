#' The Fit class
#'
#' @export
StanModelFit <- R6::R6Class("StanModelFit",
  private = list(
    stan_fit = NULL,
    stan_data = NULL,
    datasets = NULL,
    model = NULL,
    loo = NULL
  ),
  public = list(


    #' @description
    #' Get original data used to fit the model.
    #'
    #' @param dataname Name of dataset.
    get_data = function(dataname = "LON") {
      dn <- self$datanames()
      checkmate::assert_choice(dataname, dn)
      private$datasets[[dataname]]
    },

    #' @description
    #' Get names of stored data.
    #'
    datanames = function() {
      names(private$datasets)
    },

    #' @description
    #' Get model.
    #'
    #' @param name Name of model. Currently as no effect.
    get_model = function(name) {
      private$model
    },

    #' @description
    #' Create model fit object
    #'
    #' @param model A \code{\link{StanModel}} object.
    #' @param stan_fit The created 'Stan' fit.
    #' @param datasets Original data sets.
    #' @param stan_data The created 'Stan' data. Stored mainly for easier
    #' debugging.
    initialize = function(model, stan_fit, datasets, stan_data) {
      checkmate::assert_list(datasets)
      checkmate::assert_class(model, "StanModel")
      private$model <- model
      private$datasets <- datasets
      private$stan_fit <- stan_fit
      private$stan_data <- stan_data
    },

    #' @description Get the underlying 'Stan' fit object.
    get_stan_fit = function() {
      private$stan_fit
    },

    #' @description Get the underlying 'Stan' data object (for debugging).
    get_stan_data = function() {
      private$stan_data
    },

    #' @description Extract draws as \code{rvar}s
    #'
    #' @param name Param/quantity name
    draws = function(name = NULL) {
      d <- self$get_stan_fit()$draws(name)
      d <- posterior::as_draws_rvars(d)
      if (is.null(name)) {
        return(d)
      }
      d[[name]]
    },

    #' @description Get total number of draws.
    num_draws = function() {
      posterior::ndraws(self$loglik())
    },

    #' @description
    #' Extract log likelihood as 'rvars'.
    loglik = function() {
      self$draws("log_lik")
    },

    #' @description
    #' Extract the log likelihood matrix with shape \code{n_draws} x
    #' \code{n_eval_points}.
    loglik_matrix = function() {
      sf <- self$get_stan_fit()
      posterior::as_draws_matrix(sf$draws("log_lik"))
    },

    #' @description
    #' Returns estimate and standard error for ELPD
    #'
    #' @param rerun Re-run loo even if its result is already cached.
    #' @param ... Arguments passed to \code{loo:loo()}. Will only matter
    #' if re-running or running loo for the first time.
    loo_estimate = function(rerun = FALSE, ...) {
      loo <- self$loo_object(rerun = rerun, ...)
      as.numeric(loo$estimates[1, ])
    },

    #' @description
    #' Returns loo object.
    #'
    #' @param rerun Re-run loo even if its result is already cached.
    #' @param ... Arguments passed to \code{loo:loo()}. Will only matter
    #' if re-running or running loo for the first time.
    loo_object = function(rerun = FALSE, ...) {
      if (is.null(private$loo) || rerun) {
        message("loo not run yet or rerun=TRUE, running loo now...")
        log_lik <- self$loglik()
        private$loo <- loo::loo(posterior::as_draws_matrix(log_lik), ...)
      }
      private$loo
    },

    #' @description Generate quantities using the fit.
    #'
    #' @param stan_data Full 'Stan' input list.
    #' @param fitted_params Argument to \code{generate_quantities()}.
    #' @param ... Other arguments to \code{generate_quantities()}.
    gq = function(stan_data = NULL, fitted_params = NULL, ...) {
      if (is.null(fitted_params)) {
        fitted_params <- self$get_stan_fit()
      }
      if (inherits(fitted_params, "draws_rvars")) {
        fitted_params <- posterior::as_draws_array(fitted_params)
      }

      # Call 'Stan'
      self$get_model()$get_stanmodel()$generate_quantities(
        fitted_params = fitted_params,
        data = stan_data,
        ...
      )
    }
  )
)
