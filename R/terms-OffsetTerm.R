# Used by FormulaParser
create_offsetterm <- function(covariates, hierarchy) {
  z_name <- covariates[1]
  h_name <- hierarchy
  if (is.null(h_name)) {
    out <- GroupedOffsetTerm$new(z_name)
  } else {
    stop("h_name should be NULL for OffsetTerm")
  }
  out
}


# Define the OffsetTerm class
# Abstract class which should not be instantiated
OffsetTerm <- R6::R6Class("OffsetTerm",
  inherit = FormulaTerm,
  lock_class = TRUE,
  public = list(
    prior_intercept = "student_t(4, 0, 5)",
    initialize = function() {
      private$suffix <- "0"
      private$typename <- "baseline"
    },
    param_name = function() {
      ts <- private$suffix
      paste0("offset_", ts)
    },
    stancode_model = function() {
      paste0(self$param_name(), " ~ ", self$prior_intercept, ";\n")
    }
  )
)

# Define the grouped OffsetTerm class
GroupedOffsetTerm <- R6::R6Class("GroupedOffsetTerm",
  inherit = OffsetTerm,
  lock_class = TRUE,
  public = list(
    initialize = function(z_name) {
      super$initialize()
      checkmate::assert_character(z_name)
      self$z_name <- z_name
      private$suffix <- z_name
    },
    stancode_data = function(used_names, datanames) {
      z <- self$stanlines_data_z(datanames)
      build_stancode_data(z$lines, z$stannames, used_names)
    },
    stancode_pars = function() {
      pn <- self$param_name()
      gn <- paste0("G_", self$z_name)
      paste0("  vector[", gn, "] ", pn, ";\n")
    },
    stancode_tpars = function(datanames) {
      tn <- self$stanname(datanames)
      pn <- self$param_name()
      dnz <- self$stanname_z(datanames)
      code <- paste0(
        "  vector[n_", datanames, "] ", tn, " = ", pn, "[", dnz, "];"
      )
      paste(code, collapse = "\n")
    },
    standata = function(datasets, conf) {
      conf <- self$ensure_conf(conf)
      self$create_standata_z(datasets)
    },
    as_gam_term = function(L, B) {
      self$z_name
    }
  )
)
