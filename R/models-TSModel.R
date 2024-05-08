#' Can be used to just generate Stan code for a TS model
# given a formula
#'
#' @export
#' @param formula model formula
#' @param print Should the code be printed?
#' @param ... Arguments passed to the \code{$create_stancode()} method
#' of \code{\link{TSModel}}.
#' @return Stan code as character string (invisibly)
stancode_ts <- function(formula, print = TRUE, ...) {
  checkmate::assert_class(formula, "formula")
  m <- TSModel$new(formula, compile = F)
  code <- m$create_stancode(...)
  if (print) {
    cat(code)
  }
  invisible(code)
}

#' Time series (longitudinal) model class (R6 class)
#'
#' @export
#' @field term_list The additive model terms.
#' @field y_var Name of the y variable.
#' @field id_var Name of the subject identifier variable.
#' @field prior_sigma Prior for the noise parameter.
#' @field sigma_upper Upper bound for the noise parameter.
#' @field sigma_lower Lower for the noise parameter.
#' @field log_y_cap Upper bound on log scale for creating a capped
#' predicted signal or predicted observations.
TSModel <- R6::R6Class("TSModel",
  inherit = StanModel,

  # PRIVATE
  private = list(
    loglik_suffix = "lon",
    default_dataname = function(datanames) {
      if (is.null(datanames)) {
        datanames <- "LON"
      }
      datanames
    },
    stanfiles_functions_impl = function() {
      self$term_list$stanfiles_functions()
    },
    stancode_data_impl = function(datanames) {
      datanames <- private$default_dataname(datanames)
      code <- paste0(
        "  int<lower=1> n_", datanames, ";",
        " // number of input points (", datanames, ")"
      )
      paste(
        paste(code, collapse = "\n"),
        self$term_list$stancode_data(datanames),
        stancode_ts_data(self$stanname_y(), "LON"),
        sep = "\n"
      )
    },
    stancode_tdata_impl = function(datanames) {
      datanames <- private$default_dataname(datanames)
      dn_def <- private$default_dataname(NULL)
      c1 <- stancode_ts_tdata(self$stanname_y(), dn_def)
      c2 <- self$term_list$stancode_tdata(datanames)
      paste(c1, c2, sep = "\n")
    },
    stancode_pars_impl = function() {
      scb_sigma <- stancode_bounds(self$sigma_lower, upper = self$sigma_upper)
      c1 <- paste0("  real", scb_sigma, " sigma; // noise magnitude \n")
      c2 <- self$term_list$stancode_pars()
      paste(c1, c2, sep = "\n")
    },
    stancode_tpars_impl = function(datanames) {
      datanames <- private$default_dataname(datanames)
      code <- self$term_list$stancode_tpars(datanames)
      dn_def <- private$default_dataname(NULL)
      paste(code,
        stancode_ts_likelihood(self$stanname_y(), dn_def),
        sep = "\n"
      )
    },
    stancode_model_impl = function() {
      c1 <- paste0("  sigma ~ ", self$prior_sigma, ";\n")
      c2 <- self$term_list$stancode_model()
      c3 <- stancode_loglik(private$loglik_suffix)
      paste(c1, c2, c3, sep = "\n")
    },
    stancode_gq_impl = function() {
      terms_gq <- self$term_list$stancode_gq()
      stancode_ts_gq(self, self$stanname_y(), terms_gq)
    },
    delta = NULL
  ),

  # PUBLIC
  public = list(
    id_var = NULL,
    term_list = NULL,
    y_var = NULL,
    prior_sigma = "normal(0, 2)",
    sigma_upper = 3,
    sigma_lower = 0,
    log_y_cap = 7,

    #' @description
    #' Create model
    #'
    #' @param formula The model formula determining the terms and the y
    #' variable (longitudinal observation).
    #' @param id_var Name of the subject identifier variable.
    #' @param compile Should the 'Stan' model code be created and compiled.
    #' @param delta Offset for log transform (\code{y_log = log(y + delta)}).
    #' @param prior_sigma Prior for sigma
    #' @param prior_terms A list with names equal to a subset of the
    #' names of the model terms. Can be used to edit priors of term parameters.
    #' @param prior_baseline Prior for the baseline term.
    #' @param sigma_upper Upper bound for sigma
    #' @param sigma_lower Lower bound for sigma
    #' @param baseline Baseline term definition. Created automatically based
    #' on \code{id_var} if \code{NULL} (default).
    initialize = function(formula, id_var = "id", compile = TRUE, delta = 0,
                          baseline = NULL,
                          prior_baseline = NULL,
                          prior_terms = NULL,
                          prior_sigma = "normal(0, 2)",
                          sigma_upper = 3,
                          sigma_lower = 0) {
      checkmate::assert_character(id_var, min.chars = 1)
      checkmate::assert_class(formula, "formula")

      # Handle adding the baseline term
      complete <- complete_formula_ts(formula, id_var, baseline)
      formula <- complete$formula
      baseline <- complete$baseline
      prior_terms <- complete_prior_terms(prior_terms, prior_baseline, baseline)

      # Set fields
      self$prior_sigma <- prior_sigma
      self$sigma_upper <- sigma_upper
      self$sigma_lower <- sigma_lower
      self$term_list <- create_termlist(formula, prior_terms)
      self$y_var <- FormulaParser$new(formula)$get_y_name()
      self$id_var <- id_var
      private$delta <- delta

      # Compile
      super$initialize(compile)
    },

    #' @description
    #' Get name of y variable in Stan code.
    stanname_y = function() {
      paste0("dat_", self$y_var)
    },

    #' @description
    #' The model description as a string
    string = function() {
      str <- paste0(
        class_name_hl(self), ":\n    ",
        "log(", self$y_var, "+", self$get_delta(), ") ~ N(f, sigma^2)"
      )
      tls <- self$term_list$string()
      paste0(str, ", where \n", "    f = ", hl_string(tls))
    },

    #' @description Get value of \code{delta}.
    get_delta = function() {
      private$delta
    },

    #' @description Get term names in Stan code.
    #' @return A character vector with length equal to number of terms.
    term_names = function() {
      self$term_list$stan_names()
    },

    #' @description
    #' Create the 'Stan' data list from a data frame. Performs normalization
    #' on continuous variables that are input to GP terms.
    #'
    #' @param data A data frame.
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @param num_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param scale_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param skip_transform Term names whose input transform should be
    #' skipped.
    #' @param prior_only Sample from prior only?
    #' @param set_transforms If data transforms should be set based on the given
    #' \code{data}. This should be \code{TRUE} when fitting a model, and
    #' \code{FALSE} when computing predictions using GQ.
    #' @param dataname Name of dataset.
    #' @return A list.
    create_standata = function(data,
                               term_confs = NULL,
                               num_bf = NULL,
                               scale_bf = NULL,
                               skip_transform = NULL,
                               prior_only = FALSE,
                               set_transforms = TRUE,
                               dataname = "LON") {
      checkmate::assert_data_frame(data, col.names = "named")
      checkmate::assert_logical(set_transforms)
      checkmate::assert_character(dataname)

      # Fill missing term configurations
      full_term_confs <- self$term_list$fill_term_confs(
        term_confs,
        num_bf = num_bf,
        scale_bf = scale_bf
      )

      # Set variable transforms
      if (set_transforms) {
        self$term_list$set_transforms(data, skip_transform)
      }

      # Prepare 'Stan' data
      data <- ensure_id_var_exists(data, self$id_var)

      # Create final data list
      stan_data <- standata_ts(
        data, dataname, self$term_list, full_term_confs, self$y_var,
        private$delta
      )

      # Return
      list(
        orig_data = data,
        stan_data = finalize_stan_data(stan_data, prior_only),
        full_term_confs = full_term_confs
      )
    },


    #' @description
    #' Fit the model.
    #'
    #' @param data A data frame.
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @param num_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param scale_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param skip_transform Term names whose input transform should be
    #' skipped.
    #' @param prior_only Sample from prior only.
    #' @param ... Arguments passed to \code{sample} method of the
    #' 'CmdStanR' model.
    #' @return An \code{\link{TSModelFit}} object.
    fit = function(data,
                   term_confs = NULL,
                   num_bf = NULL,
                   scale_bf = NULL,
                   skip_transform = NULL,
                   prior_only = FALSE,
                   ...) {
      # Get Stan model object
      stan_model <- self$get_stanmodel()
      if (is.null(stan_model)) {
        stop("Stan model does not exist, you need to call compile()!")
      }

      # Create Stan input list
      d <- self$create_standata(
        data, term_confs, num_bf, scale_bf,
        skip_transform, prior_only, TRUE
      )

      # Call 'Stan'
      stan_fit <- stan_model$sample(data = d$stan_data, ...)

      # Return
      dat_list <- list(LON = d$orig_data)
      TSModelFit$new(self, stan_fit, dat_list, d$stan_data, d$full_term_confs)
    }
  )
)


# Stan data for TS model
standata_ts <- function(data, dataname, term_list, term_confs, y_name, delta) {
  sd_y <- standata_ts_y(data, y_name, delta)
  sd_x <- term_list$create_standata(data, dataname, term_confs)
  sd <- c(sd_x, sd_y)
  sd$prior_only <- 0
  sd
}

# Stan data specific to TS model
standata_ts_y <- function(data, y_name, delta) {
  checkmate::assert_character(y_name)
  y <- data[[y_name]]
  checkmate::assert_numeric(y)
  out <- list(y, delta)
  names(out) <- c(paste0("dat_", y_name), "delta")
  out
}

# Ensure that data has id specifier
ensure_id_var_exists <- function(data, id_var) {
  id_var_given <- id_var %in% colnames(data)
  msg <- paste0(id_var, " not found in data, assuming data are from same id")
  if (!(id_var_given)) {
    message(msg)
    data[[id_var]] <- as.factor(rep(1, nrow(data)))
  }
  data
}

# Complete the formula for TS model
complete_formula_ts <- function(formula, id_var, baseline) {
  ff <- as.character(formula)
  if (is.null(baseline)) {
    baseline <- paste0("offset(", id_var, ")")
    message("baseline was NULL, setting baseline = ", hl_string(baseline))
  }
  checkmate::assert_character(baseline, min.chars = 1)
  if (ff[3] == ".") {
    formula <- as.formula(paste0(ff[2], " ~ ", baseline))
  } else {
    dpf <- deparse(formula)
    dpf <- paste(dpf, collapse = "")
    formula <- as.formula(paste(dpf, "+", baseline))
  }
  list(
    formula = formula,
    baseline = baseline
  )
}

# Complete the prior list for terms, handle baseline separately
complete_prior_terms <- function(prior_terms, prior_baseline, baseline) {
  baseline_stan <- term_to_code(baseline)
  if (baseline_stan %in% names(prior_terms)) {
    stop(
      "prior for baseline should not be given in prior_terms, instead use ",
      "prior_baseline"
    )
  }
  prior_terms[[baseline_stan]] <- prior_baseline
  prior_terms
}
