# Base class
ForwardSearch <- R6Class(
  classname = "ForwardSearch",
  public = list(
    num_comps = NULL,
    num_steps = NULL,
    verbosity = 0,
    score_name = "unknown",
    path = NULL, # possible predefined path
    draw_inds = NULL,
    draw_inds_eval = NULL,
    # Init
    initialize = function(J, draw_inds, draw_inds_eval, path = NULL,
                          num_steps = 1) {
      checkmate::assert_integerish(path, null.ok = TRUE)
      checkmate::assert_integerish(num_steps, lower = 1)
      checkmate::assert_integerish(J, lower = 1)
      checkmate::assert_integerish(draw_inds, lower = 1)
      checkmate::assert_integerish(draw_inds_eval, lower = 1)
      self$num_comps <- J
      self$verbosity <- 1
      self$path <- path
      self$num_steps <- num_steps
      self$draw_inds <- draw_inds
      self$draw_inds_eval <- draw_inds_eval
    },

    # Score a model
    score = function(model, ...) {
      stop("Inheriting class should override the score method!\n")
    },

    # Compute score for each candidate
    step = function(model, candidates, ...) {
      stop("Inheriting class should override the step method!\n")
    },

    # Get candidate models at current step
    get_candidates = function(cur_model, J) {
      j <- length(cur_model) + 1
      if (length(self$path) >= j) {
        cands <- self$path[j] # using pre-defined path
        message("path has been defined for this step")
      } else {
        cands <- setdiff(seq_len(J), cur_model)
        message("path not defined for this step -> possibly many candidates")
      }
      cands
    }
  )
)

# Forward search using projection
ProjectionForwardSearch <- R6Class(
  classname = "ProjectionForwardSearch",
  inherit = ForwardSearch,
  public = list(
    score_name = "negative_KL",

    # Score a model
    score = function(model, fit_ref, eval_mode, kl0, elpd_loo_ref) {
      checkmate::assert_integerish(model, unique = TRUE, null.ok = TRUE)
      if (!eval_mode) {
        di <- self$draw_inds
      } else {
        di <- self$draw_inds_eval
      }
      res <- fit_ref$project(model, eval_mode = eval_mode, draw_inds = di)
      res <- res$metrics
      res$score <- -res$kl
      if (is.null(kl0)) {
        kl0 <- res$kl
      }
      res$p_explained <- 1.0 - res$kl / kl0
      edr <- as.numeric((elpd_loo_ref[1] - res$elpd_loo) / elpd_loo_ref[2])
      res$elpd_loo_rel_diff <- edr
      res
    },

    # Choose best candidate
    step = function(model, candidates, fit_ref) {
      if (length(candidates) == 1) {
        return(1)
      }
      df <- NULL
      elpd_loo_ref <- fit_ref$loo_estimate()
      for (ca in candidates) {
        c_model <- c(model, ca)
        res <- self$score(c_model, fit_ref, FALSE, NULL, elpd_loo_ref)
        df <- rbind(df, res)
      }
      df <- data.frame(df)
      cand_scores <- df$score
      which(cand_scores == max(cand_scores))
    },

    # Perform search
    run = function(fit_ref, ...) {
      J <- self$num_comps
      model <- NULL
      score <- NULL
      elpd <- NULL
      elpd_loo_ref <- fit_ref$loo_estimate()
      history <- self$score(model, fit_ref, TRUE, NULL, elpd_loo_ref)
      kl0 <- history$kl
      j <- 0

      # Loop
      while (length(model) < J) {
        j <- j + 1
        msg <- paste0("Step ", j, "/", J, ".")
        message(msg)

        # Choose best candidate
        cands <- self$get_candidates(model, J)
        message("candidate terms: {", paste0(cands, collapse = ", "), "}")
        i_best <- self$step(model, cands, fit_ref)

        # Score the best candidate
        new_model <- c(model, cands[i_best])
        new_row <- self$score(new_model, fit_ref, TRUE, kl0, elpd_loo_ref)
        print(new_row)

        # Update history
        history <- rbind(history, new_row)
        model <- c(model, cands[i_best])
        edr <- new_row$elpd_loo_rel_diff
        cat("elpd_loo_rel_diff = ", edr, "\n", sep = "")
        if (is.null(self$num_steps)) {
          if (edr < 1) {
            break
          }
        } else {
          if (j >= self$num_steps) {
            break
          }
        }
      }

      # Return
      list(path = model, history = history)
    }
  )
)


#' Run projection predictive forward search
#'
#' @export
#' @param fit_ref The reference model fit.
#' @param path Pre-defined search path. If \code{NULL}, path is
#' not pre-defined. These indices are same as the index of the term in
#' \code{fit_ref$get_model()$term_names()}.
#' @param num_steps Number of steps to run. If \code{NULL}, the search stops
#' when predictive performance of the submodel is close to that of reference
#' model.
#' @param S Number of draws to project when selecting next step.
#' @param S_eval Number of draws to project when evaluating step.
pp_forward_search <- function(fit_ref, path = NULL, num_steps = NULL,
                              S = 30, S_eval = 100) {
  checkmate::assert_class(fit_ref, "LonModelFit")
  J <- length(fit_ref$get_model()$term_list$terms)
  S_tot <- fit_ref$num_draws()
  message(
    "Using ", S, " (", S_eval, ") out of ", S_tot, " draws for ",
    "search (evaluation)"
  )
  draw_inds <- sample.int(S_tot, S)
  draw_inds_eval <- sample.int(S_tot, S_eval)
  a <- ProjectionForwardSearch$new(
    J, draw_inds, draw_inds_eval, path, num_steps
  )
  a$run(fit_ref)
}

#' Plot result of projection predictive forward search (p_explained)
#'
#' @export
#' @param res The list returned by \code{\link{pp_forward_search}}
#' @param thresh A horizontal line.
plot_pp_pexp <- function(res, thresh = 0.95) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  p_exp <- res$history$p_explained
  df <- data.frame(num_vars, p_exp)
  out <- ggplot(df, aes(x = num_vars, y = p_exp)) +
    geom_hline(yintercept = thresh, color = "firebrick3", lty = 2) +
    ylab("1 - KL / KL0")
  out <- out + geom_line() + geom_point() +
    xlab("Number of terms") + ggtitle("Forward search")
  return(out)
}

#' Plot result of projection predictive forward search (ELPD)
#'
#' @export
#' @param res The list returned by \code{\link{pp_forward_search}}
#' @param elpd_ref ELPD of reference model (estimate and SE)
plot_pp_elpd <- function(res, elpd_ref) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  elpd <- res$history$elpd_loo
  elpd_se <- res$history$elpd_loo_se
  df <- data.frame(num_vars, elpd, elpd_se)
  out <- ggplot(df, aes(
    x = num_vars, y = elpd, ymin = elpd - elpd_se,
    ymax = elpd + elpd_se
  )) +
    geom_errorbar(width = 0.1) +
    geom_hline(yintercept = elpd_ref[1], color = "steelblue3", lty = 3) +
    geom_hline(
      yintercept = elpd_ref[1] - elpd_ref[2], color = "steelblue3",
      lty = 1
    ) +
    geom_hline(
      yintercept = elpd_ref[1] + elpd_ref[2], color = "steelblue3",
      lty = 1
    ) +
    ylab("ELPD")
  elpd <- res$history$elpd
  st <- paste0(
    "blue = ref. model PSIS-LOO, ",
    "black = submodel PSI-LOO, ",
    "red = submodel full data"
  )
  out <- out + geom_line() + geom_point() +
    geom_line(data = data.frame(num_vars, elpd), lty = 2, color = "firebrick3") +
    xlab("Number of terms") +
    ggtitle("Forward search", subtitle = st)
  return(out)
}
