# Used by FormulaParser
create_offsetterm <- function(covariates, hierarchy) {
  z_name <- covariates[1]
  h_name <- hierarchy
  if (is.null(h_name)) {
    out <- GroupedOffsetTerm$new(z_name)
  } else {
    out <- HierOffsetTerm$new(z_name = z_name, h_name = h_name)
  }
  out
}


# Define the OffsetTerm class
# Abstract class which should not be instantiated
OffsetTerm <- R6::R6Class("OffsetTerm",
  inherit = FormulaTerm,
  lock_class = TRUE,
  private = list(
    latex_type = "BAS"
  ),
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
  private = list(
    latex_param_names = function() {
      "\\mathbf{c}"
    }
  ),
  public = list(
    initialize = function(z_name) {
      super$initialize()
      checkmate::assert_character(z_name)
      self$z_name <- z_name
      private$suffix <- z_name
      private$latex_param_subscript <- "0"
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
    }
  )
)


# Define the grouped OffsetTerm class
HierOffsetTerm <- R6::R6Class("HierOffsetTerm",
  inherit = OffsetTerm,
  lock_class = TRUE,
  private = list(
    latex_param_names = function() {
      "\\mathbf{c}"
    }
  ),
  public = list(
    prior_log_z = "std_normal()",
    prior_log_mu = "std_normal()",
    prior_log_sigma = "normal(0, 5)",
    initialize = function(z_name, h_name) {
      super$initialize()
      checkmate::assert_character(z_name)
      checkmate::assert_character(h_name)
      self$z_name <- z_name
      self$h_name <- h_name
      private$suffix <- z_name
      private$latex_param_subscript <- "0"
    },
    stancode_data = function(used_names, datanames) {
      z <- self$stanlines_data_z(datanames)
      h <- self$stanlines_data_h(datanames)
      lines <- c(z$lines, h$lines)
      sn <- c(z$stannames, h$stannames)
      build_stancode_data(lines, sn, used_names)
    },
    stancode_pars = function() {
      gz <- paste0("G_", self$z_name)
      gh <- paste0("G_", self$h_name)
      raw_params_stancode(self$param_name(), gh, gz)
    },
    stancode_tpars = function(datanames) {
      zn <- self$stanname_z(datanames)
      hn <- self$stanname_h(datanames)
      pn <- self$param_name()
      fc <- stancall_offset_grouped(pn, hn, zn)
      sn <- self$stanname(datanames)
      code <- paste0(
        "  vector[n_", datanames, "] ", sn, " = ", fc, ";"
      )
      paste(code, collapse = "\n")
    },
    stancode_model = function() {
      prior <- list(
        log_z = self$prior_log_z,
        log_mu = self$prior_log_mu,
        log_sigma = self$prior_log_sigma
      )
      hier_prior_stancode(self$param_name(), prior)
    },
    standata = function(datasets, conf) {
      conf <- self$ensure_conf(conf)
      sd_z <- self$create_standata_z(datasets)
      sd_h <- self$create_standata_h(datasets)
      c(sd_z, sd_h)
    }
  )
)

# Parameter to natural scale (in transformed parameters)
stancall_offset_grouped <- function(parname, h_name, z_name) {
  codes <- c()
  J <- length(h_name)
  for (j in seq_len(J)) {
    ofs <- param_vecs_log_hier(parname, h_name[j], z_name[j])
    in_par <- paste(ofs, collapse = ", ")
    code_j <- paste0("to_log_natural_scale(", in_par, ")")
    codes <- c(codes, code_j)
  }
  codes
}
