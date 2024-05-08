#' Create additive terms
#'
#' @param formula the model formula
#' @param priors Prior configurations.
#' @return An object of \code{R6} class \code{TermList}.
create_termlist <- function(formula, priors) {
  checkmate::assert_class(formula, "formula")
  fp <- FormulaParser$new(formula)
  terms <- fp$parse_terms()
  list <- TermList$new(terms)
  pr_names <- names(priors)
  for (pr in pr_names) {
    t <- list$get_term(pr)
    message("Editing term ", hl_string(pr))
    prior <- priors[[pr]]
    checkmate::assert_list(prior)
    for (nam in names(prior)) {
      message(" * Editing field ", hl_string(nam))
      t[[nam]] <- prior[[nam]]
    }
  }
  list
}

#' Additive terms (R6 class)
#'
#' @export
#' @field terms A list of \code{FormulaTerm} objects, defining the
#' function components.
#' @field fsum_name Base name of the sum of the terms variable in 'Stan' code.
#' @description A semi-parametric additive model (R6 class).
TermList <- R6::R6Class("TermList",
  inherit = StanCodeCreator,

  # PRIVATE
  private = list(
    stan_model = NULL,
    finalize_code = function(code) {
      if (nchar(trimws(code)) == 0) {
        return("")
      }
      checkmate::assert_character(code, len = 1)
      paste0(
        "\n  // ", class_name(self), " for ", self$fsum_name, "\n",
        code, "\n"
      )
    },
    stanfiles_functions_impl = function() {
      gtr <- function(x) x$stanfiles_functions()
      unlist(sapply(self$terms, gtr))
    },
    stancode_data_impl = function(datanames) {
      stancode_termlist_data(self$terms, datanames)$code
    },
    stancode_tdata_impl = function(datanames) {
      f <- function(x) x$finalize_code(x$stancode_tdata(datanames))
      paste(sapply(self$terms, f), collapse = "\n")
    },
    stancode_pars_impl = function() {
      f <- function(x) x$finalize_code(x$stancode_pars())
      paste(sapply(self$terms, f), collapse = "\n")
    },
    stancode_tpars_impl = function(datanames) {
      # Terms
      f1 <- function(x) x$finalize_code(x$stancode_tpars(datanames))
      sc_tpars <- paste(sapply(self$terms, f1), collapse = "\n")

      # Sum
      sc_sum <- stancode_f_sum(datanames, self$terms, self$fsum_name)
      paste(sc_tpars, sc_sum, sep = "\n")
    },
    stancode_model_impl = function() {
      f <- function(x) x$finalize_code(x$stancode_model())
      paste(sapply(self$terms, f), collapse = "\n")
    },
    stancode_gq_impl = function() {
      dataname <- "LON"
      fsum_name <- self$fsum_name
      nn <- paste0("n_", dataname)
      f <- function(x) x$finalize_code(x$stancode_gq())
      sc_gq <- paste(sapply(self$terms, f), collapse = "\n")
      paste0(
        sc_gq, "\n  // Sum of terms\n  vector[", nn, "] ", fsum_name,
        " = ", fsum_name, "_", dataname, ";\n"
      )
    }
  ),

  # PUBLIC
  public = list(
    terms = NULL,
    fsum_name = "f_sum",

    #' @description
    #' Create term list
    #'
    #' @param terms The list of model terms.
    initialize = function(terms) {
      checkmate::assert_class(terms, "list")
      self$terms <- terms
    },

    #' @description
    #' The term list as a string
    string = function() {
      gn <- function(x) x$stanname_base()
      paste0(sapply(self$terms, gn), collapse = " + ")
    },

    #' @description
    #' Get names of all input variables
    input_vars = function() {
      getter <- function(x) x$input_vars()
      inputs_list <- sapply(self$terms, getter)
      vars <- unique(unlist(inputs_list))
    },

    #' @description List names of all terms in Stan code.
    stan_names = function() {
      gtr <- function(x) x$stanname_base()
      out <- sapply(self$terms, gtr)
      names(out) <- NULL
      out
    },

    #' @description Get term object based on term name in Stan'
    #' @param f_name_stan Name of term in 'Stan'. Valid names
    #' can be checked using the `$stan_names()` method.
    get_term = function(f_name_stan) {
      sn <- self$stan_names()
      idx <- find_term_index(f_name_stan, sn, TRUE)
      if (length(idx) != 1) {
        stop("Internal error. Please report a bug.")
      }
      self$terms[[idx]]
    },

    #' @description
    #' Latex code, mathematical notation for the terms
    latex = function() {
      J <- length(self$terms)
      fs <- paste0("f^{(", seq_len(J), ")}(\\mathbf{x})")
      math <- sapply(self$terms, function(x) x$latex())
      latex_terms <- paste0(fs, " &= ", math)
      paste(
        "\\begin{align}\n  ",
        paste(latex_terms, collapse = "\\\\ \n  "),
        "\\end{align}",
        sep = "\n"
      )
    },

    #' @description List length
    length = function() {
      length(self$terms)
    },

    #' @description
    #' Create the 'Stan' data list from a data frame. Performs transforms
    #' on continuous variables that are input to GP terms.
    #' @param data A data frame.
    #' @param dataname Name of data set.
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @return A list.
    create_standata = function(data, dataname, term_confs = NULL) {
      term_names <- names(self$terms)
      sd <- list(nrow(data))
      names(sd) <- paste0("n_", dataname) # dimension
      datasets <- list(data)
      names(datasets) <- dataname
      for (tn in term_names) {
        sd_term <- self$terms[[tn]]$standata(datasets, term_confs[[tn]])
        sd <- c(sd, sd_term)
      }
      sd[unique(names(sd))] # remove duplicates
    },

    #' @description
    #' Get term configuration of all terms, some of which can be given by user.
    #'
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @param num_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param scale_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @return The updated model object (invisibly).
    fill_term_confs = function(term_confs = NULL,
                               num_bf = NULL, scale_bf = NULL) {
      stannames <- names(self$terms)
      confs <- list()
      j <- 0
      for (tn in stannames) {
        j <- j + 1
        term <- self$terms[[tn]]
        tc <- term$ensure_conf(term_confs[[tn]])
        if (is(term, "FormulaSFTerm")) {
          tc$kg <- term$term_list_kg$fill_term_confs(tc, num_bf, scale_bf)
          tc$ks <- term$term_list_ks$fill_term_confs(tc, num_bf, scale_bf)
        }
        if (is(term, "GPTerm")) {
          tc <- fill_term_conf_gpterm(tc, num_bf, scale_bf)
        }
        confs[[j]] <- tc
      }
      names(confs) <- stannames
      confs
    },

    #' @description
    #' Set variable transforms using data.
    #'
    #' @param data A data frame.
    #' @param names_skip Names of terms whose input transform should be
    #' skipped.
    set_transforms = function(data, names_skip = NULL) {
      settable <- c("Normalizing", "MaxScale", "UnitScale")
      settable <- paste0(settable, "Transform")
      ts <- self$terms
      for (tn in names(ts)) {
        term <- ts[[tn]]
        transf <- term$x_transform
        msg <- paste0("Setting transform for ", hl_string(tn), "... ")
        if (tn %in% names_skip) {
          msg <- paste0(msg, "skipped!")
        } else {
          if (inherits(transf, settable)) {
            x_data <- get_x_from_data(data, term$x_name)
            term$x_transform <- transf$set_using_data(x_data)
            transf_name <- class(transf)[1]
            msg <- paste0(msg, "done (", hl_string(transf_name), ")!")
          } else {
            msg <- paste0(msg, "not needed.")
          }
        }
        message(msg)

        # Handle FormulaSFTerm
        if (inherits(term, "FormulaSFTerm")) {
          message("Setting transforms for the terms inside of FormulaSFTerm")
          term$term_list_kg$set_transforms(data, names_skip)
          term$term_list_ks$set_transforms(data, names_skip)
        }
      }
    },


    #' @description
    #' Check how data transforms using the current model transforms.
    #'
    #' @param data A data frame.
    check_transforms = function(data) {
      settable <- c("Normalizing", "MaxScale", "UnitScale")
      settable <- paste0(settable, "Transform")
      ts <- self$terms
      for (tn in names(ts)) {
        term <- ts[[tn]]
        transf <- term$x_transform
        msg <- paste0("Checking transform for ", hl_string(tn), "... ")
        if (inherits(transf, settable)) {
          x_data <- get_x_from_data(data, term$x_name)
          r1 <- range(x_data)
          r2 <- range(transf$forward(x_data))
          transf_name <- class(transf)[1]
          msg <- paste0(msg, "(", hl_string(transf_name), ")!")
          message(msg)
          message("Data range is [", r1[1], ", ", r1[2], "]")
          message("Transformed range is [", r2[1], ", ", r2[2], "]")
        } else {
          msg <- paste0(msg, "no transform.")
          message(msg)
        }
      }
    }
  )
)

# Term configuration of GPTerm
fill_term_conf_gpterm <- function(tc, num_bf, scale_bf) {
  if (!is.null(num_bf)) {
    tc$num_bf <- num_bf
  }
  if (!is.null(scale_bf)) {
    tc$scale_bf <- scale_bf
  }
  tc
}

# Find index of model term based on Stan variable name
find_term_index <- function(f_name_stan, term_names_stan, strict) {
  sn <- term_names_stan
  tn_str <- paste(sn, collapse = ", ")
  if (!(f_name_stan %in% sn) && strict) {
    msg <- paste0(
      hl_string(f_name_stan), " is not a Stan variable name of any",
      " term in this model. Valid names are {", hl_string(tn_str), "}."
    )
    stop(msg)
  }
  which(sn == f_name_stan)
}


# Generate data block
stancode_termlist_data <- function(terms, datanames, code = NULL,
                                   used_data_names = NULL) {
  for (t in terms) {
    scd <- t$stancode_data(used_names = used_data_names, datanames)
    used_data_names <- unique(c(used_data_names, scd$stannames))
    code_add <- t$finalize_code(scd$code)
    code <- paste(code, code_add, sep = "\n")
  }
  list(code = code, used_names = used_data_names)
}

# Helper
stancode_f_sum <- function(datanames, terms, fsum_name) {
  fun <- function(x) x$stanname(datanames)
  f_list <- lapply(terms, fun) # gives list with length = num terms
  J <- length(datanames)
  code <- "  // Sum the terms"
  for (j in seq_len(J)) {
    nn <- paste0("n_", datanames[j])
    gtr <- function(x) x[j]
    f_names <- paste(lapply(f_list, gtr), collapse = " + ")
    code_j <- paste0(
      "  vector[", nn, "] ", fsum_name, "_", datanames[j], " = ", f_names, ";"
    )
    code <- paste(code, code_j, sep = "\n")
  }
  code
}


# Get sff term of a model
get_sff_term <- function(model) {
  terms <- model$term_list$terms
  cl <- sapply(terms, function(x) class(x)[1] == "FormulaSFTerm")
  idx <- which(cl)
  if (length(idx) != 1) {
    stop("The model should have one FormulaSFTerm")
  }
  terms[[which(cl)]]
}
