# Define the GPTerm class
GPTerm <- R6::R6Class("GPTerm",
  inherit = FormulaTerm,
  lock_class = TRUE,
  private = list(
    stanfiles_functions_impl = function() {
      c("gp/eigenfunctions", "gp/eigenvalues")
    }
  ),
  public = list(
    x_transform = UnitScaleTransform$new(),
    initialize = function(x_name, z_name = NULL) {
      if (is.na(z_name)) {
        z_name <- NULL
      }
      checkmate::assert_character(x_name, any.missing = FALSE)
      if (!is.null(z_name)) {
        checkmate::assert_character(z_name, any.missing = FALSE)
      }
      self$x_name <- x_name
      self$z_name <- z_name
      private$typename <- "gp"
      private$suffix <- stancode_suffix_interact(x_name, z_name)
    },
    default_conf = function() {
      list(scale_bf = 1.5, num_bf = 30)
    },
    required_conf_names = function() {
      c("scale_bf", "num_bf")
    },
    stancode_data = function(used_names, datanames) {
      x <- self$stanlines_data_x(datanames)
      gp <- stanlines_data_hsgp(private$suffix)
      if (self$has_z()) {
        z <- self$stanlines_data_z(datanames)
      } else {
        z <- list(lines = NULL, stannames = NULL)
      }
      lines <- c(x$lines, gp$lines, z$lines)
      names <- c(x$stannames, gp$stannames, z$stannames)
      build_stancode_data(lines, names, used_names)
    },
    stancode_tdata = function(datanames) {
      tn <- private$suffix
      xn <- self$stanname_x(datanames)
      bn <- paste0("B_", tn)
      mn <- paste0("mat_B_", tn, "_", datanames)
      sn <- paste0("seq_B_", tn)
      c1 <- paste0("  vector[", bn, "] ", sn, " = seq_len(", bn, ");\n")
      c2 <- paste0("  ", stancode_B_matrix(datanames, bn, mn, sn), ";\n")
      if (self$has_z()) {
        zn <- self$stanname_z(datanames)
        gn <- paste0("G_", self$z_name)
        c3 <- paste0("  ", stancode_PSI_grouped(datanames, xn, tn, bn, mn, gn, zn), ";\n")
      } else {
        c3 <- paste0("  ", stancode_PSI_shared(datanames, xn, tn, bn, mn), ";\n")
      }

      paste0(c(c1, c2, c3), collapse = "")
    },
    stancode_pars = function() {
      tn <- private$suffix
      gn <- paste0("G_", self$z_name)
      code <- paste0("\n  real<lower=1e-12> alpha_", tn, "; // magnitude")
      code <- paste0(code, "\n  real<lower=1e-12> ell_", tn, "; // lengthscale")
      if (!self$has_z()) {
        code <- paste0(code, "\n  vector[B_", tn, "] xi_", tn, "; // auxiliary")
      } else {
        code <- paste0(code, "\n  vector[(", gn, " - 1) * B_", tn, "] xi_", tn, ";")
      }
      paste0(code, "\n")
    },
    stancode_tpars = function(datanames) {
      sfx <- private$suffix
      sn <- self$stanname(datanames)
      mn <- paste0("seq_B_", sfx)
      bn <- paste0("B_", sfx)
      if (self$has_z()) {
        gn <- paste0("G_", self$z_name)
        c1 <- stancall_DELTA_grouped(sfx, bn, mn, gn)
      } else {
        c1 <- stancall_DELTA_shared(sfx, bn, mn)
      }
      c2 <- stancall_gp_eval(datanames, sn, sfx)
      paste0(c(c1, c2), collapse = "")
    },
    stancode_model = function() {
      tn <- private$suffix
      code <- paste0("\n  alpha_", tn, " ~ student_t(20, 0, 1);")
      code <- paste0(code, "\n  ell_", tn, " ~ lognormal(0, 1);")
      code <- paste0(code, "\n  xi_", tn, " ~ normal(0, 1);")
      paste0(code, "\n")
    },
    standata = function(datasets, conf) {
      datanames <- names(datasets)
      conf <- self$ensure_conf(conf)
      ts <- private$suffix

      # Stan data for continuous variable
      xn <- self$x_name
      dat <- self$create_standata_x(datasets)
      sn_x <- self$stanname_x(datanames)

      # HSGP Stan data
      term_names <- self$stanname(datanames)
      j <- 0
      for (tn in term_names) {
        j <- j + 1
        dat_x_unit <- dat[[sn_x[j]]]
        dat_hsgp_j <- create_standata_hsgp(ts, conf, dat_x_unit, xn, tn)
        dat <- c(dat, dat_hsgp_j)
      }

      # Stan data for categorical variable
      if (self$has_z()) {
        dat_z <- self$create_standata_z(datasets)
      } else {
        dat_z <- NULL
      }

      # Return
      c(dat, dat_z)
    },
    as_gam_term = function(L, B) {
      s <- create_hs_smooth(self$stanname_x("LON"), self$stanname_z("LON"), L, B)
      z <- self$stanname_z("LON")
      if (!self$has_z()) {
        return(s)
      }
      paste0(s, " + ", z)
    }
  )
)

# Helper
stancode_B_matrix <- function(datasets, bn, mn, sn) {
  nn <- paste0("n_", datasets)
  paste0(
    "matrix[", nn, ", ", bn, "] ", mn, " = transpose(rep_matrix(",
    sn, ", ", nn, "))"
  )
}


# Helper
stancode_PSI_shared <- function(datasets, xn, tn, bn, mn) {
  nn <- paste0("n_", datasets)
  tn_sfx <- paste0(tn, "_", datasets)
  paste0(
    "matrix[", nn, ", ", bn, "] PSI_", tn_sfx,
    " = eigfun_shared(", xn, ", ", mn, ", L_", tn, ")"
  )
}

# Helper
stancode_PSI_grouped <- function(datasets, xn, tn, bn, mn, gn, zn) {
  nn <- paste0("n_", datasets)
  tn_sfx <- paste0(tn, "_", datasets)
  paste0(
    "matrix[", nn, ", ", bn, " * (", gn, " - 1)] PSI_", tn_sfx,
    " = eigfun_grouped(", xn, ", ", mn, ", L_", tn,
    ", ", zn, ", ", gn, " )"
  )
}


# Helper
stancall_DELTA_shared <- function(tn, bn, sn) {
  paste0(
    "vector[", bn, "] DELTA_", tn,
    " = eigval_shared(alpha_", tn, ", ell_", tn, " , ", sn, ", L_", tn, ");\n"
  )
}

# Helper
stancall_DELTA_grouped <- function(tn, bn, sn, gn) {
  paste0(
    "vector[", bn, " * (", gn, " - 1)] DELTA_", tn,
    " = eigval_grouped(alpha_", tn, ", ell_", tn, " , ", sn, ", L_", tn,
    ", ", gn, " );\n"
  )
}

# Prevent improper use
ensure_gp_approx_is_valid <- function(dat_x_unit, x_name, L, term_name) {
  max_abs_x <- max(abs(dat_x_unit))
  x_name <- paste0(x_name, " (transformed to unit scale)")
  if (max_abs_x > L) {
    msg <- paste0("Error in term '", term_name, "':\n")
    msg <- paste0(
      msg,
      " max absolute value of ", x_name, " is larger than L = ", L,
      ".\n GP approximation not valid for this input."
    )
    stop(msg)
  }
}

# Stan data for HSGP approximation
create_standata_hsgp <- function(suffix, conf, dat_x_unit, x_name, term_name) {
  checkmate::assert_number(conf$scale_bf, lower = 0)
  checkmate::assert_integerish(conf$num_bf, lower = 1)
  L <- conf$scale_bf
  ensure_gp_approx_is_valid(dat_x_unit, x_name, L, term_name)
  out <- list(L, conf$num_bf)
  names(out) <- paste0(c("L_", "B_"), suffix)
  out
}

# Helper
stanlines_data_hsgp <- function(suffix) {
  dn <- list(
    L = paste0("L_", suffix),
    B = paste0("B_", suffix)
  )
  lines <- list(
    paste0("  real<lower=0> ", dn$L, "; // domain size"),
    paste0("  int<lower=1> ", dn$B, "; // num of basis funs")
  )
  list(
    lines = lines,
    stannames = unlist(dn)
  )
}

# Shared GP term
stancall_gp_eval <- function(datanames, sn, sfx) {
  nn <- paste0("n_", datanames)
  sfx_full <- paste0(sfx, "_", datanames)
  paste0(
    "  vector[", nn, "] ", sn,
    " = PSI_", sfx_full, " * (xi_", sfx, " .* DELTA_", sfx, ");\n"
  )
}

# Group-specific GP term
stancall_gp_group <- function(datanames, sn, sfx, zn) {
  nn <- paste0("n_", datanames)
  sfx_full <- paste0(sfx, "_", datanames)
  paste0(
    "  vector[", nn, "] ", sn,
    " = compute_f_group(xi_", sfx, ", PHI_", sfx_full, ", s_", sfx,
    ", VARPHI_", sfx_full, ", d_", sfx, ");\n"
  )
}
