#' The 'LonModelFit' class
#'
#' @export
#' @field term_confs The term configurations used.
LonModelFit <- R6::R6Class("LonModelFit",
  inherit = StanModelFit,
  private = list(

    # Extract one component as 'rvar'
    extract_f_comp = function(f_name = "f") {
      self$draws(name = f_name)
    }
  ),
  public = list(
    term_confs = NULL,

    #' @description
    #' Create model fit object
    #'
    #' @param model A \code{\link{StanModel}} object.
    #' @param stan_fit The created 'Stan' fit.
    #' @param datasets Original data (list of data frames).
    #' @param stan_data The created 'Stan' data. Stored mainly for easier
    #' debugging.
    #' @param term_confs The term configurations used.
    initialize = function(model, stan_fit, datasets, stan_data, term_confs) {
      self$term_confs <- term_confs
      super$initialize(model, stan_fit, datasets, stan_data)
    },

    #' @description Print object description.
    print = function() {
      cat("An R6 object of class LonModelFit.\n")
    },

    #' @description
    #' Extract one component as a \code{FunctionDraws} object.
    #'
    #' @param component Name of the component (term). Valid names can
    #' be checked using \code{$model$term_names()}. Alternatively, can be
    #' an integer that is the component index.
    #' @param dataname name of data set used to evaluate the function
    function_draws = function(component = "f_sum") {
      mod <- self$get_model("lon")
      dat <- self$get_data("LON")
      if (is.numeric(component)) {
        checkmate::assert_integerish(component)
        component <- mod$term_names()[component]
      }
      f_name <- component

      if (f_name == "f_sum" || f_name == "y_log_pred") {
        if (f_name == "y_log_pred") {
          f_name <- paste0(mod$y_var, "_log_pred")
        }
        # The total sum or predictive
        f <- private$extract_f_comp(f_name)
        covs <- mod$term_list$input_vars()
        x <- dat[, covs, drop = FALSE]
        fd <- FunctionDraws$new(x, f, f_name)
      } else {
        # A single component
        t <- mod$term_list$get_term(f_name)
        covs <- t$input_vars()
        x <- dat[, covs, drop = FALSE]
        f <- private$extract_f_comp(f_name)
        fd <- FunctionDraws$new(x, f, f_name)
      }
      fd
    },

    #' @description
    #' Initialize a \code{\link{FunctionDraws}} object from the fit.
    #'
    #' @param input_vars Column names in the data.
    #' @param f_name Name of function in 'Stan' code.
    create_functiondraws = function(input_vars, f_name) {
      dat <- self$get_data("LON")
      x <- dat[, input_vars]
      f <- private$extract_f_comp(f_name)
      FunctionDraws$new(x, f, f_name)
    },

    #' @description
    #' Plot fit or predictive distribution
    #'
    #' @param f_reference Reference function to plot against the fit.
    #' @param plot_y Should y data be plotted?
    #' @param f_ref_name Name of the reference function.
    #' @param predictive Should this plot the predictive distribution for
    #' new observations? Otherwise the inferred signal is plotted.
    #' @param filter_by Factor to filter the rows by.
    #' @param kept_vals Values of the factor that are not filtered out.
    #' @param ... Arguments passed to the plot method of
    #' \code{\link{FunctionDraws}}.
    plot = function(f_reference = NULL, plot_y = TRUE,
                    predictive = TRUE,
                    f_ref_name = "the true signal",
                    filter_by = NULL,
                    kept_vals = NULL,
                    ...) {
      dat <- self$get_data("LON")
      yn <- self$get_model("lon")$y_var
      if (predictive) {
        f_name <- "y_log_pred"
      } else {
        f_name <- "f_sum"
      }
      f <- self$function_draws(component = f_name)
      plt <- f$plot(filter_by = filter_by, kept_vals = kept_vals, ...)
      x_var <- rlang::as_name(plt$mapping$x)
      dat <- df_filter_rows(dat, filter_by, kept_vals)

      # Plot true signal as comparison
      if (!is.null(f_reference)) {
        checkmate::assert_numeric(f_reference, len = nrow(dat))
        dat$f_reference <- f_reference
        message("black dashed line is ", f_ref_name)
        plt <- plt + geom_line(aes(x = !!sym(x_var), y = f_reference),
          color = "black", lty = 2,
          data = dat,
          inherit.aes = FALSE
        )
      }

      # Plot data as dots
      if (plot_y) {
        plt <- plt + geom_point(
          aes(x = !!sym(x_var), y = !!sym(yn)),
          color = "black",
          data = dat,
          inherit.aes = FALSE
        )
      }
      plt
    },

    #' Compute predictions at test points
    #'
    #' @param data The input data frame of test points.
    #' @param fitted_params Parameters or model fit.
    #' @param ... Other arguments to \code{generate_quantities()}.
    #' @return A new fit object.
    predict = function(data = NULL, fitted_params = NULL, ...) {
      # Create input
      if (is.null(data)) {
        data <- self$get_data("LON")
      }
      model <- self$get_model()

      # Create Stan input list
      dat <- model$create_standata(
        data = data,
        term_confs = self$term_confs,
        num_bf = NULL,
        scale_bf = NULL,
        skip_transform = NULL,
        prior_only = FALSE,
        set_transforms = FALSE
      )

      # Call 'Stan'
      gq <- self$gq(
        stan_data = dat$stan_data, fitted_params = fitted_params,
        ...
      )

      # Return
      data_orig <- list(LON = dat$orig_data)
      LonModelFit$new(
        model, gq, data_orig, dat$stan_data, dat$full_term_confs
      )
    },

    #' Compute predictions at given time points
    #'
    #' @param t_test The test time points
    #' @param t_var Name of the time variable
    #' @param t_var_copy Names of other continuous covariates to which
    #' \code{t_test} should be copied.
    #' @param y_test Test y
    #' @param ... Arguments passed to \code{$predict()}.
    #' @return A new fit object
    predict_time = function(t_test, t_var = "t", y_test = NULL,
                            t_var_copy = NULL, ...) {
      dat_lon <- self$get_data("LON")
      y_name <- self$get_model()$y_var
      df_test <- create_df_predict_time(
        dat_lon, t_test, t_var, t_var_copy,
        y_test, y_name
      )
      self$predict(df_test, ...)
    }
  )
)

# Helper
create_df_predict_time <- function(dat_lon, t_test, t_var, t_var_copy,
                                   y_test, y_name) {
  df_test <- extend_df(dat_lon, t_test, t_var)
  for (v in t_var_copy) {
    df_test[[v]] <- df_test[[t_var]]
  }
  if (is.null(y_test)) {
    y_test <- 1
  }
  df_test[[y_name]] <- y_test
  df_test
}