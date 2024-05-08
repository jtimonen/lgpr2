#' Class for objects that can create parts of 'Stan' code
#'
StanCodeCreator <- R6::R6Class("StanCodeCreator",
  private = list(
    finalize_code = function(code) {
      if (nchar(trimws(code)) == 0) {
        return("")
      }
      checkmate::assert_character(code, len = 1)
      paste0("\n  // ", class_name(self), "\n", code, "\n")
    },
    stanfiles_functions_impl = function() {
      NULL
    },
    stancode_data_impl = function(datanames) {
      ""
    },
    stancode_tdata_impl = function(datanames) {
      ""
    },
    stancode_pars_impl = function() {
      ""
    },
    stancode_tpars_impl = function(datanames) {
      ""
    },
    stancode_model_impl = function() {
      ""
    },
    stancode_gq_impl = function() {
      ""
    }
  ),
  public = list(

    #' @description
    #' Description as a string.
    string = function() {
      paste0(class_name_hl(self))
    },

    #' @description
    #' Print info.
    #'
    #' @return Nothing.
    print = function() {
      cat(self$string(), "\n", sep = "")
    },

    #' @description Stan files where to get the functions.
    stanfiles_functions = function() {
      files <- c("utils", private$stanfiles_functions_impl())
      unique(files)
    },

    #' @description Generate 'Stan' code for the functions block.
    stancode_functions = function() {
      files <- self$stanfiles_functions()
      code <- ""
      for (fn in files) {
        code <- paste0(code, read_stan_function(fn), sep = "\n")
      }
      code
    },

    #' @description Generate 'Stan' code for the data block.
    #' @param datanames Names of input data sets.
    stancode_data = function(datanames = NULL) {
      private$finalize_code(private$stancode_data_impl(datanames))
    },

    #' @description Generate 'Stan' code for the transformed data block.
    #' @param datanames Names of input data sets.
    stancode_tdata = function(datanames = NULL) {
      private$finalize_code(private$stancode_tdata_impl(datanames))
    },

    #' @description Generate 'Stan' code for the parameters block.
    stancode_pars = function() {
      private$finalize_code(private$stancode_pars_impl())
    },

    #' @description Generate 'Stan' code for the transformed parameters block.
    #' @param datanames Names of input data sets.
    stancode_tpars = function(datanames = NULL) {
      private$finalize_code(private$stancode_tpars_impl(datanames))
    },

    #' @description Generate 'Stan' code for the model block.
    stancode_model = function() {
      private$finalize_code(private$stancode_model_impl())
    },

    #' @description Generate 'Stan' code for the generated quantities block.
    stancode_gq = function() {
      private$finalize_code(private$stancode_gq_impl())
    }
  )
)


#' Abstract for objects that can generate entire 'Stan' models
#'
StanModel <- R6::R6Class("StanModel",
  inherit = StanCodeCreator,

  # PRIVATE
  private = list(
    stan_model = NULL,
    stancode_tpars_impl = function(datanames) {
      paste0("  real log_lik_empty = 0.0;")
    }
  ),

  # PUBLIC
  public = list(

    #' @description
    #' Create model
    #'
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(compile = TRUE) {
      if (compile) {
        self$compile()
      }
    },

    #' @description Get the underlying 'Stan' model.
    get_stanmodel = function() {
      private$stan_model
    },

    #' @description
    #' Create the 'Stan' model code for the model.
    #'
    #' @param autoformat Should automatic formatting be attempted?
    #' @return A string.
    create_stancode = function(autoformat = TRUE) {
      checkmate::check_logical(autoformat)
      sc_data <- paste(
        "  int<lower=0,upper=1> prior_only;\n",
        self$stancode_data()
      )
      code <- paste(
        stancode_block("functions", self$stancode_functions()),
        stancode_block("data", sc_data),
        stancode_block("transformed data", self$stancode_tdata()),
        stancode_block("parameters", self$stancode_pars()),
        stancode_block("transformed parameters", self$stancode_tpars()),
        stancode_block("model", self$stancode_model()),
        stancode_block("generated quantities", self$stancode_gq()),
        sep = "\n"
      )

      if (autoformat) {
        tryCatch(
          {
            code <- autoformat_stancode(code)
          },
          error = function(e) {
            cat(code)
            stop(e)
          }
        )
      }
      code
    },

    #' @description
    #' Create and compile the 'Stan' model.
    #'
    #' @param dir Path to directory where to store the \code{.stan} file.
    #' @return A 'CmdStanR' model.
    create_stanmodel = function(dir = tempdir()) {
      code <- self$create_stancode(autoformat = FALSE)
      a <- cmdstanr::write_stan_file(code = code, dir = dir)

      # silence compile warnings from cmdstan
      utils::capture.output(
        {
          mod <- cmdstanr::cmdstan_model(a)
        },
        type = "message"
      )
      mod
    },

    #' @description
    #' Create and compile the 'Stan' model.
    #'
    #' @param ... Arguments passed to \code{create_stanmodel()}.
    #' @return The updated model object (invisibly).
    compile = function(...) {
      private$stan_model <- self$create_stanmodel(...)
      invisible(self)
    }
  )
)

# Helper
stancode_block <- function(name, code) {
  paste0(name, "{\n", code, "\n}\n")
}
