# Define the purely abstract base class FormulaTerm
FormulaTerm <- R6::R6Class("FormulaTerm",
  inherit = StanCodeCreator,
  lock_class = TRUE,
  private = list(
    typename = NULL,
    suffix = NULL
  ),
  public = list(
    x_name = NULL,
    z_name = NULL,
    h_name = NULL, # hierarchy
    x_transform = IdentityTransform$new(),
    has_x = function() {
      !is.null(self$x_name)
    },
    has_z = function() {
      !is.null(self$z_name)
    },
    has_hierarchy = function() {
      !is.null(self$h_name)
    },
    input_vars = function() {
      unique(c(self$x_name, self$z_name, self$h_name))
    },
    finalize_code = function(code) {
      checkmate::assert_character(code, len = 1)
      code <- trimws(code)
      if (nchar(code) == 0) {
        return("")
      }
      tn_long <- self$name_long()
      paste0("  // Term ", tn_long, "\n  ", code, "\n")
    },
    stanname_x = function(datanames) {
      x_name <- paste0("dat_", self$x_name)
      paste0(self$x_transform$add_suffix(x_name), "_", datanames)
    },
    stanname_z = function(datanames) {
      paste0("dat_", self$z_name, "_", datanames)
    },
    stanname_h = function(datanames) {
      paste0("dat_", self$h_name, "_", datanames)
    },
    stanname_base = function() {
      paste0("f_", private$typename, "_", private$suffix)
    },
    stanname = function(datanames) {
      paste0(self$stanname_base(), "_", datanames)
    },
    name_long = function() {
      paste0(self$stanname_base(), " (", class_name(self), ")")
    },
    stancode_data = function(used_names, datanames) {
      list(code = "", stannames = character(0))
    },
    stancode_gq = function() {
      dataname <- "LON"
      tn <- self$stanname_base()
      tn_tpar <- self$stanname(dataname)
      nn <- paste0("n_", dataname)
      paste0("  vector[", nn, "] ", tn, " = ", tn_tpar, ";")
    },
    print = function() {
      cat(class_name(self), "(base name in Stan code = ",
        self$stanname_base(), ").\n",
        sep = ""
      )
    },
    as_gam_term = function(L, B) {
      stop("as_gam_term() not implemented for this term type")
    },
    standata = function(datasets, conf) {
      list()
    },
    stanlines_data_x = function(datanames) {
      dnx <- self$stanname_x(datanames)
      nn <- paste0("n_", datanames)
      lines <- as.list(
        paste0("  vector[", nn, "] ", dnx, "; // continuous input")
      )
      list(lines = lines, stannames = c(dnx))
    },
    stanlines_data_z = function(datanames) {
      dnz <- self$stanname_z(datanames)
      dnG <- paste0("G_", self$z_name)
      nn <- paste0("n_", datanames)
      zlines <- paste0("  array[", nn, "] int ", dnz, "; // categ input")
      lines <- c(
        as.list(zlines),
        list(paste0("  int<lower=1> ", dnG, "; // num of groups"))
      )
      list(lines = lines, stannames = c(dnz, dnG))
    },
    create_standata_x = function(datasets) {
      checkmate::assert_list(datasets, min.len = 1)
      datanames <- names(datasets)
      j <- 0
      OUT <- list()
      for (df in datasets) {
        j <- j + 1
        checkmate::assert_data_frame(df, col.names = "named")
        x <- get_x_from_data(df, self$x_name)
        x <- self$x_transform$forward(x) # always give transformed data to Stan
        out <- list(x)
        snx <- self$stanname_x(datanames[j])
        names(out) <- c(snx)
        OUT <- c(OUT, out)
      }
      OUT
    },
    create_standata_z = function(datasets) {
      checkmate::assert_list(datasets, min.len = 1)
      nn <- names(datasets)
      create_standata_categorical(datasets, self$z_name, self$stanname_z(nn))
    },
    default_conf = function() {
      list()
    },
    required_conf_names = function() {
      c()
    },

    # this should not be overridden by inheriting class
    ensure_conf = function(conf) {
      if (length(conf) == 0) {
        conf <- self$default_conf()
      }
      checkmate::assert_list(conf, names = "named")
      cn <- self$required_conf_names()
      for (nam in cn) {
        if (!(nam %in% names(conf))) {
          stop(nam, " not found in conf list!")
        }
      }
      conf
    }
  )
)

# Util
name_does_not_exist <- function(name, names = NULL) {
  checkmate::assert_character(name)
  if (!is.null(names)) {
    checkmate::assert_character(names)
  }
  !(name %in% names)
}

# Appends line to Stan code if variable name is not already in use
append_to_datablock <- function(code, line, var_name, used_names) {
  if (name_does_not_exist(var_name, used_names)) {
    code <- paste0(code, line, "\n")
  }
  code
}

# To be called in the end of stancode_data methods of FormulaTerms
build_stancode_data <- function(lines, stannames, used_names) {
  code <- ""
  J <- length(stannames)
  if (length(lines) != J) {
    stop(
      "length(lines) != length(stannames) in build_stancode_data(), ",
      "please report a bug."
    )
  }
  for (j in seq_len(J)) {
    code <- append_to_datablock(code, lines[[j]], stannames[j], used_names)
    used_names <- c(used_names, stannames[j])
  }
  list(
    code = code,
    stannames = stannames
  )
}

# Get continuous variable named 'xn' from data frame
get_x_from_data <- function(data, xn) {
  x <- data[[xn]]
  if (is.null(x)) {
    stop(xn, " not found in data!")
  }
  checkmate::assert_numeric(x)
  x
}

# Helper
create_standata_categorical <- function(datasets, name, stannames) {
  zn <- name
  j <- 0
  OUT <- list()
  for (df in datasets) {
    j <- j + 1
    checkmate::assert_data_frame(df, col.names = "named")
    z <- df[[zn]]
    if (is.null(z)) {
      stop(zn, " not found in data!")
    }
    checkmate::assert_factor(z)
    G <- length(levels(z))
    out <- list(
      as.numeric(z),
      G
    )
    gn <- paste0("G_", zn)
    nams <- c(stannames[j], gn)
    names(out) <- nams
    OUT <- c(OUT, out)
  }
  OUT
}

# Stan code variable name suffix for interaction term
stancode_suffix_interact <- function(x_name, z_name) {
  ts <- x_name
  if (!is.null(z_name)) {
    ts <- paste0(ts, "X", z_name)
  }
  ts
}
