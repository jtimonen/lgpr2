#' The 'LonModelFit' class
#'
#' @export
#' @field term_confs The term configurations used.
LonModelFit <- R6::R6Class("LonModelFit",
  inherit = StanModelFit,
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

    #' @description Stan variable names.
    #' @param all Return all Stan variable names. If \code{FALSE}, returns only
    #' ones that are not large transformed parameters.
    stan_variables = function(all = TRUE) {
      md <- self$get_stan_fit()$metadata()
      sv <- md$stan_variables
      if (all) {
        return(sv)
      }
      pattern <- "^(?!f_|DELTA_).*$"
      grep(pattern, sv, value = TRUE, perl = TRUE)
    },

    #' @description
    #' A diagnostic summary.
    diagnose = function() {
      sf <- self$get_stan_fit()
      vars <- self$stan_variables(all = FALSE)
      rhat <- sf$summary(vars)$rhat
      diag <- unlist(sf$diagnostic_summary())
      out <- c(diag, max(rhat))
      names(out)[length(out)] <- "max_rhat"
      out
    },

    #' @description
    #' Extract one component as a \code{FunctionDraws} object.
    #'
    #' @param component Name of the component (term). Valid names can
    #' be checked using \code{$model$term_names()}. Alternatively, can be
    #' an integer that is the component index.
    #' @param input_vars Column names in the data. Only has effect if
    #' getting a single component. If \code{NULL}, determined automatically
    #' from the component.
    #' @param dataname name of data set used to evaluate the function
    function_draws = function(component = "h", input_vars = NULL) {
      mod <- self$get_model("lon")
      dat <- self$get_data("LON")
      if (is.numeric(component)) {
        checkmate::assert_integerish(component)
        component <- mod$term_names()[component]
      }
      f_name <- component

      if (f_name == "f_sum" || f_name == "y_log_pred" || f_name == "h") {
        if (f_name == "y_log_pred") {
          f_name <- paste0(mod$y_var, "_log_pred")
        }

        # The total sum or predictive
        f <- self$draws(f_name)
        covs <- mod$term_list$input_vars()

        x <- dat[, covs, drop = FALSE]
        fd <- FunctionDraws$new(x, f, f_name)
      } else {
        # A single component
        t <- mod$term_list$get_term(f_name)
        if (is.null(input_vars)) {
          covs <- t$input_vars()
        } else {
          covs <- input_vars
        }
        x <- dat[, covs, drop = FALSE]
        f <- self$draws(f_name)
        fd <- FunctionDraws$new(x, f, f_name)
      }
      fd
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
        f_name <- "h"
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
        set_transforms = FALSE,
        set_c_hat = FALSE
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
        dat_lon, t_test, t_var, t_var_copy, y_test, y_name
      )
      self$predict(df_test, ...)
    },

    #' Estimate amount of noise
    #'
    #' @return A tibble with num of rows equal to number of draws. Values
    #' are between 0 (no noise) and 1 (data is only noise).
    noise_amount = function() {
      m <- self$get_model()
      y_data <- self$get_data()[[m$y_var]]
      fd <- self$function_draws()
      h <- fd$as_data_frame_long()
      h %>%
        dplyr::group_by(.draw_idx) %>%
        dplyr::summarize(p_noise = compute_noise_amount(y_data, value))
    },

    #' Estimate component variances for each draw
    #'
    #' @return A matrix with n_draws rows and n_comps columns
    comp_vars = function() {
      m <- self$get_model()
      cn <- m$term_names()
      out <- list()
      for (comp_idx in 1:length(cn)) {
        fd <- self$function_draws(component = cn[comp_idx])
        df <- fd$as_data_frame_long() %>% dplyr::group_by(.draw_idx)
        df <- df %>% dplyr::summarize(var = stats::var(value))
        out[[comp_idx]] <- df$var
      }
      names(out) <- cn
      sapply(out, function(x) x)
    },

    #' Estimate component relevances
    #'
    #' @param rvar return an \code{rvar} vector?
    relevances = function(rvar = TRUE) {
      p_noise <- self$noise_amount()$p_noise
      p_comp <- self$comp_vars()
      p_comp <- p_comp / rowSums(p_comp) * (1 - p_noise)
      out <- cbind(p_comp, p_noise)
      if (rvar) {
        out <- posterior::rvar(out)
      }
      out
    },


    #' @description
    #' Rank model terms based on relevance
    #'
    rank_terms = function() {
      p_exp <- mean(self$relevances())
      J <- length(p_exp)
      p_noise <- p_exp[J]
      p_exp <- p_exp[1:(J - 1)]
      srt <- sort(p_exp, decreasing = TRUE, index.return = TRUE)
      rels <- c(p_noise, srt$x)
      names(rels)[1] <- "p_noise"
      list(
        relevances = rels,
        order = srt$ix
      )
    },


    #' Reduce model
    #'
    #' @param thresh Threshold for explained variance.
    reduce = function(thresh = 0.95) {
      p_exp <- self$rank_terms()
      num_sel <- length(which(cumsum(p_exp$relevances) < thresh)) + 1
      rels <- p_exp$relevances
      list(
        relevance = p_exp$relevances,
        order = p_exp$order,
        num_sel = num_sel,
        selected = names(rels)[seq_len(num_sel)]
      )
    },

    #' Project to submodel
    #'
    #' @param term_inds Which terms to include in submodel?
    #' @param draw_inds Which posterior draws to use. If \code{NULL}, 30
    #' draws are taken randomly.
    #' @param eval_mode Model evaluation mode?
    #' @param B number of basis functions
    #' @param random_idx Index of random term. If used, the projection
    #' is a GAMM instead of a GAM.
    project = function(term_inds = NULL, draw_inds = NULL, eval_mode = TRUE,
                       B = 24, random_idx = NULL) {
      m <- self$get_model()
      term_inds <- setdiff(term_inds, random_idx)
      form <- m$as_gam_formula(term_inds = term_inds, B = B)
      if (!is.null(random_idx)) {
        random_form <- m$term_list$as_random_offset_formula(random_idx)
      } else {
        random_form <- NULL
      }
      print(form)
      print(random_form)


      # Prepare
      h_ref <- self$function_draws()
      S <- h_ref$num_draws()
      if (is.null(draw_inds)) {
        if (eval_mode) {
          S <- 100
        } else {
          S <- 30
        }
        draw_inds <- sample.int(S, size = S)
      }
      h_df <- h_ref$as_data_frame_long()
      h_df <- h_df %>% dplyr::filter(.draw_idx %in% draw_inds)
      dat <- self$get_data()
      dat <- transform_df(m, dat)

      # Do the work
      pd <- project_draws(self, dat, h_df, form, random_form)

      # Metrics
      loglik_mat <- pd$loglik_proj
      elpd <- mean(colSums(loglik_mat))
      if (eval_mode) {
        loo <- loo::loo(t(loglik_mat))
        loo_est <- as.numeric(loo$estimates[1, ])
      } else {
        loo <- NULL
        loo_est <- c(NA, NA)
      }

      # Return
      metrics <- data.frame(
        kl = mean(pd$kl_div),
        elpd = elpd,
        elpd_loo = loo_est[1],
        elpd_loo_se = loo_est[2]
      )
      list(
        metrics = metrics,
        mu_proj = pd$mu_proj,
        mu_ref = pd$mu_ref,
        num_fails = 0,
        loglik = t(loglik_mat),
        loo = loo
      )
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
