#' The 'FunctionDraws' class
#'
#' @export
FunctionDraws <- R6::R6Class("FunctionDraws",
  private = list(
    input = NULL,
    output = NULL,
    name = NULL,
    var_types = function() {
      sapply(private$input, class)
    },
    cont_vars = function() {
      idx_cont <- which(private$var_types() == "numeric")
      if (length(idx_cont) == 0) {
        return(NULL)
      }
      colnames(private$input)[idx_cont]
    },
    categ_vars = function() {
      idx_cat <- which(private$var_types() == "factor")
      if (length(idx_cat) == 0) {
        return(NULL)
      }
      colnames(private$input)[idx_cat]
    },

    # Create data frame for ggplot
    quantiles_df_full = function(ci_inner, ci_outer) {
      get_q <- function(x, p) {
        f <- as.vector(stats::quantile(x, probs = p))
      }
      p_in_low <- (1 - ci_inner) / 2
      p_in_up <- 1 - p_in_low
      p_out_low <- (1 - ci_outer) / 2
      p_out_up <- 1 - p_out_low
      f_draws <- private$output
      out_low <- get_q(f_draws, p_out_low)
      in_low <- get_q(f_draws, p_in_low)
      med <- get_q(f_draws, 0.5)
      in_up <- get_q(f_draws, p_in_up)
      out_up <- get_q(f_draws, p_out_up)
      f_dist <- cbind(out_low, in_low, med, in_up, out_up)
      cbind(private$input, data.frame(f_dist))
    }
  ),
  public = list(

    #' @description Create the object.
    #' @param input The data frame used as function input.
    #' @param output An \code{rvar} containing draws of the function values
    #' at \code{input} points.
    #' @param name Function name.
    initialize = function(input, output, name) {
      checkmate::assert_character(name, min.chars = 1)
      checkmate::assert_class(input, "data.frame")
      checkmate::assert_class(output, "rvar")
      checkmate::assert_true(length(dim(input)) == 2)
      checkmate::assert_true(length(dim(output)) == 1)
      checkmate::assert_true(nrow(input) == length(output))
      private$input <- input
      private$output <- output # rvar
      private$name <- name
    },

    #' @description
    #' Get number of function draws
    num_draws = function() {
      posterior::ndraws(private$output)
    },

    #' @description
    #' Get number of input points.
    num_points = function() {
      length(private$output)
    },

    #' @description
    #' Get the function name.
    get_name = function() {
      private$name
    },

    #' @description
    #' Rename the function.
    #'
    #' @param name new name
    rename = function(name) {
      checkmate::assert_character(name, len = 1)
      private$name <- name
    },

    #' @description
    #' Get the input data frame.
    get_input = function() {
      private$input
    },

    #' @description
    #' Get the output \code{rvar}.
    #'
    #' @param as_matrix Transform the \code{rvar} to a draws matrix
    #' (see the \code{posterior} package).
    get_output = function(as_matrix = FALSE) {
      out <- private$output
      if (as_matrix) {
        out <- posterior::as_draws_matrix(out)
      }
      out
    },

    #' @description
    #' Get the input and output as \code{data.frame} with an \code{rvar}
    #' column.
    as_data_frame = function() {
      value <- self$get_output()
      cbind(self$get_input(), value)
    },

    #' @description
    #' Get the input and output as a long \code{data.frame} with each draw
    #' stacked on top of each other.
    as_data_frame_long = function() {
      S <- self$num_draws()
      N <- self$num_points()
      row_idx <- rep(1:N, each = S)
      draw_idx <- rep(1:S, times = N)
      df <- self$get_input()[row_idx, ]
      df[[".draw_idx"]] <- as.factor(draw_idx)
      df$value <- as.vector(self$get_output(as_matrix = TRUE))
      df
    },

    #' @description
    #' Create (a possibly filtered) data frame for ggplot.
    #' @param ci_inner inner credible interval
    #' @param ci_outer outer credible interval
    #' @param filter_by filtering variable
    #' @param kept_vals value of the filter variable that are kept
    quantiles_df = function(ci_inner = 0.5, ci_outer = 0.8,
                            filter_by = NULL, kept_vals = NULL) {
      df <- private$quantiles_df_full(ci_inner = ci_inner, ci_outer = ci_outer)
      df_filter_rows(df, filter_by, kept_vals)
    },

    #' @description Print object description.
    print = function() {
      S <- self$num_draws()
      P <- self$num_points()
      cat("An R6 object of class FunctionDraws (",
        number_string(S), " draws, ", number_string(P), " points). \n",
        sep = ""
      )
      v1 <- private$cont_vars()
      s1 <- paste(v1, collapse = ", ")
      cat(" - Name:", private$name, "\n")
      if (length(v1) > 0) {
        cat(" - Continuous inputs:", hl_string(s1), "\n")
      }
      v2 <- private$categ_vars()
      s2 <- paste(v2, collapse = ", ")
      if (length(v2) > 0) {
        cat(" - Categorical inputs:", hl_string(s2), "\n")
      }
    },


    #' @description
    #' Plot using automatically selected aesthetics by default.
    #'
    #' @param x_var Name of the x-axis variable.
    #' @param group_by Name of the factor used as group aesthetic.
    #' @param color_by Name of the factor used as color and fill aesthetics.
    #' @param ribbon_alpha Opacity of ribbon.
    #' @param facet_by Faceting factor.
    #' @param ci_inner Inner credible interval (ribbon).
    #' @param ci_outer Outer credible interval (ribbon).
    #' @param max_levels Maximum number of levels to facet by.
    #' @param filter_by Factor to filter the rows by.
    #' @param kept_vals Values of the factor that are not filtered out.
    plot = function(x_var = "auto",
                    group_by = "auto",
                    color_by = "auto",
                    facet_by = "auto",
                    ribbon_alpha = 0.3,
                    ci_inner = 0.5,
                    ci_outer = 0.8,
                    max_levels = 40,
                    filter_by = NULL,
                    kept_vals = NULL) {
      checkmate::assert_number(ci_inner, lower = 0, upper = 1)
      checkmate::assert_number(ci_outer, lower = 0, upper = 1)
      checkmate::assert_true(ci_outer >= ci_outer)
      checkmate::assert_number(ribbon_alpha, lower = 0, upper = 1)
      df <- self$quantiles_df(
        ci_inner = ci_inner, ci_outer = ci_outer,
        filter_by = filter_by, kept_vals = kept_vals
      )
      df$Function <- as.factor(rep(private$name, nrow(df)))
      z <- private$categ_vars()[1]
      is_disabled <- function(x) {
        if (is.null(x) || is.na(x) || base::isFALSE(x)) {
          return(TRUE)
        }
        FALSE
      }
      is_auto <- function(x) {
        if (is_disabled(x)) {
          return(FALSE)
        }
        x == "auto"
      }

      # Handle auto
      if (is_auto(x_var)) {
        x_var <- private$cont_vars()[1]
      }
      if (is_auto(group_by)) {
        if (is.null(z)) {
          group_by <- NULL
        } else {
          group_by <- z
        }
      }
      if (is_auto(color_by)) {
        if (is.null(z)) {
          color_by <- NULL
        } else {
          color_by <- z
        }
      }
      if (is_auto(facet_by)) {
        if (is.null(z)) {
          facet_by <- NULL
        } else {
          facet_by <- z
        }
      }

      # Avoid slow plots
      num_z_facet <- length(unique(df[, facet_by]))
      num_z_color <- length(unique(df[, color_by]))
      max_levs_plot <- max(c(num_z_color, num_z_facet))
      if (max_levs_plot > max_levels) {
        message(
          "more than ", max_levels,
          " (max_levels) levels to color or facet by,",
          " disabling color, faceting, and ribbon"
        )
        color_by <- NULL
        facet_by <- NULL
        ribbon_alpha <- 0
      }

      # Dummy factor if NULL
      if (is_disabled(group_by)) {
        group_by <- "Function"
      }
      if (is_disabled(color_by)) {
        color_by <- "Function"
      }

      # Mapping
      mapping <- aes(
        x = !!sym(x_var), y = med,
        group = !!sym(group_by), color = !!sym(color_by),
        fill = !!sym(color_by)
      )

      # Create plot
      plt <- ggplot(df, mapping)
      if (ribbon_alpha > 0) {
        plt <- plt + geom_ribbon(
          alpha = ribbon_alpha, mapping = aes(ymin = out_low, ymax = out_up)
        ) +
          geom_ribbon(
            alpha = ribbon_alpha, mapping = aes(ymin = in_low, ymax = in_up)
          )
      }
      st <- paste0("Median")
      if (ribbon_alpha > 0) {
        ci <- paste0("{", 100 * ci_inner, ", ", 100 * ci_outer, "}")
        st <- paste0(st, " and ", ci, "% Credible Intervals")
      }
      plt <- plt + geom_line() + ggtitle(private$name, subtitle = st) + ylab("")
      if (!is_disabled(facet_by)) {
        plt <- plt + facet_wrap(stats::as.formula(paste0(".~", facet_by)),
          labeller = label_both
        )
      }
      plt
    },

    #' @description
    #' Scale the draws as \code{f_new = f*a + b}
    #' @param a multiplier
    #' @param b added constant
    scale = function(a = 1, b = 0) {
      checkmate::assert_number(a)
      checkmate::assert_number(b)
      f <- private$output * a + b # arithmetic on rvar
      new_name <- paste0(private$name, "_scaled")
      FunctionDraws$new(private$input, f, new_name)
    },

    #' @description
    #' Add to another \code{FunctionDraws}. You can also use \code{f1 + f2}.
    #' @param f A \code{FunctionDraws} object.
    add = function(f) {
      checkmate::assert_class(f, "FunctionDraws")
      x <- df_union(private$input, f$get_input())
      f_sum <- private$output + f$get_output() # sum of rvars
      nam <- paste0(private$name, " + ", f$get_name())
      FunctionDraws$new(x, f_sum, nam)
    },

    #' @description
    #' Subtract another \code{FunctionDraws}. You can also use \code{f1 - f2}.
    #' @param f A \code{FunctionDraws} object.
    subtract = function(f) {
      checkmate::assert_class(f, "FunctionDraws")
      x <- df_union(private$input, f$get_input())
      f_sub <- private$output - f$get_output() # arithmetic of rvars
      nam <- paste0(private$name, " - ", f$get_name())
      FunctionDraws$new(x, f_sub, nam)
    },

    #' @description
    #' Take log of function draws.
    log = function() {
      new_name <- paste0("log(", private$name, ")")
      FunctionDraws$new(private$input, log(private$output), new_name)
    },

    #' @description
    #' Exponentiate function draws.
    exp = function() {
      new_name <- paste0("exp(", private$name, ")")
      FunctionDraws$new(private$input, exp(private$output), new_name)
    },

    #' @description
    #' Offset every function draw so that they start from zero
    #' @return a new object of class \code{\link{FunctionDraws}}
    zero_offset = function() {
      new_out <- zero_offset_fd_rvar(self$get_output())
      FunctionDraws$new(private$input, new_out, private$name)
    },

    #' @description
    #' Split to groups
    #' @return a list of objects of class \code{\link{FunctionDraws}},
    #' one for each group/id
    #' @param id_var Name of the subject identifier factor.
    split_by_id = function(id_var) {
      split_to_subjects(self, id_var)
    },

    #' @description
    #' Offset every function draw for each subject so that they start from zero
    #' @return a new object of class \code{\link{FunctionDraws}}
    #' @param id_var Name of the subject identifier factor.
    zero_offset_by_id = function(id_var) {
      checkmate::assert_character(id_var, min.chars = 1)
      splits <- self$split_by_id(id_var)
      out <- splits[[1]]$zero_offset()$get_output()
      inp <- splits[[1]]$get_input()
      id_new <- rep(names(splits)[1], nrow(inp))
      L <- length(splits)
      if (L > 1) {
        for (j in 2:L) {
          out_j <- splits[[j]]$zero_offset()$get_output()
          inp_j <- splits[[j]]$get_input()
          id_new_j <- rep(names(splits)[j], nrow(inp_j))
          out <- c(out, out_j)
          inp <- rbind(inp, inp_j)
          id_new <- c(id_new, id_new_j)
        }
      }
      inp[[id_var]] <- as.factor(id_new)
      FunctionDraws$new(inp, out, self$get_name())
    }
  )
)

#' Add
#' @export
#' @param f1 component 1
#' @param f2 component 2
#' @return the sum \code{f1 + f2}
`+.FunctionDraws` <- function(f1, f2) {
  f1$add(f2)
}

#' Subtract
#' @export
#' @param f1 component 1
#' @param f2 component 2
#' @return the difference \code{f1 - f2}
`-.FunctionDraws` <- function(f1, f2) {
  f1$subtract(f2)
}

# Split function draws to groups by id
split_to_subjects <- function(fd, id_var) {
  inp <- fd$get_input()
  out <- fd$get_output()
  levs <- unique(inp[[id_var]])
  ret <- list()
  j <- 0
  ids_found <- c()
  for (id in levs) {
    idx_rows <- which(inp[[id_var]] == id)

    if (length(idx_rows) > 0) {
      j <- j + 1
      ids_found <- c(ids_found, id)
      input <- inp[idx_rows, , drop = FALSE]
      output <- out[idx_rows]
      name <- paste0("(", fd$get_name(), ")_", id)
      ret[[j]] <- FunctionDraws$new(input = input, output = output, name = name)
    }
  }
  names(ret) <- ids_found
  ret
}

# make every function draw start from zero
zero_offset_fd_rvar <- function(fd) {
  fd <- posterior::as_draws_array(fd) # dim n_draws x n_chains x n_points
  S <- dim(fd)[1]
  C <- dim(fd)[2]
  for (s in seq_len(S)) {
    for (j in seq_len(C)) {
      x <- as.numeric(fd[s, j, ])
      fd[s, j, ] <- x - x[1]
    }
  }
  posterior::as_draws_rvars(fd)$x
}
