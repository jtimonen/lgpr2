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
    B = NULL,
    random_idx = NULL,

    # Init
    initialize = function(J, draw_inds, draw_inds_eval, path = NULL,
                          num_steps = 1, B = 24, random_idx = NULL) {
      checkmate::assert_integerish(path, null.ok = TRUE)
      checkmate::assert_integerish(J, lower = 1)
      checkmate::assert_integerish(draw_inds, lower = 1)
      checkmate::assert_integerish(draw_inds_eval, lower = 1)
      checkmate::assert_integerish(B, lower = 1)
      self$B <- B
      self$random_idx <- random_idx
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
      res <- fit_ref$project(model,
        eval_mode = eval_mode, draw_inds = di,
        B = self$B, random_idx = self$random_idx
      )
      res <- res$metrics
      res$score <- -res$kl
      if (is.null(kl0)) {
        kl0 <- res$kl
      }
      res$p_kl <- 1.0 - res$kl / kl0
      edr_alt <- res$elpd_diff / res$elpd_diff_se
      edr <- as.numeric((elpd_loo_ref[1] - res$elpd_loo) / elpd_loo_ref[2])
      res$elpd_loo_rel_diff_alt <- edr_alt
      res$elpd_loo_rel_diff <- edr
      rownames(res) <- paste(model, collapse = "-")
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
        edr_alt <- new_row$elpd_loo_rel_diff_alt
        message("elpd_loo_rel_diff = ", round(edr, 5),
          " (alt. version = ", round(edr_alt, 5), ")",
          sep = ""
        )
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
#' @param B number of basis functions
#' @param random_idx index of the "random effect" term
pp_forward_search <- function(fit_ref, path = NULL, num_steps = NULL,
                              S = 30, S_eval = 100, B = 24, random_idx = NULL) {
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
    J, draw_inds, draw_inds_eval, path, num_steps, B, random_idx
  )
  a$run(fit_ref)
}


#' Plot result of projection predictive forward search (ELPD)
#'
#' @export
#' @param res The list returned by \code{\link{pp_forward_search}}
#' @param fit Reference model fit object.
#' @param ... Arguments passed to \code{fit$loo_estimate()}.
plot_pp_elpd <- function(res, fit, ...) {
  elpd_ref <- fit$loo_estimate(...)
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  elpd <- res$history$elpd_loo
  elpd_se <- res$history$elpd_loo_se
  df <- data.frame(num_vars, elpd, elpd_se)
  out <- ggplot(df, aes(
    x = num_vars,
    y = elpd,
    ymin = elpd - elpd_se,
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
  mlpd <- res$history$mlpd
  st <- paste0(
    "blue = ref. model PSIS-LOO, ",
    "black = submodel PSI-LOO, ",
    "red = submodel full data MLPD"
  )
  out <- out + geom_line() + geom_point() +
    geom_line(
      data = data.frame(num_vars, mlpd), lty = 2, color = "firebrick3",
      aes(x = num_vars, y = mlpd)
    ) +
    xlab("Number of terms") +
    ggtitle("ELPD", subtitle = st)
  return(out)
}


#' Plot result of projection predictive forward search (ELPD rel diff)
#'
#' @export
#' @param res The list returned by \code{\link{pp_forward_search}}
plot_pp_elpd_diff <- function(res) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  edr1 <- res$history$elpd_loo_rel_diff
  edr2 <- res$history$elpd_loo_rel_diff_alt
  df <- data.frame(num_vars = c(num_vars, num_vars), value = c(edr1, edr2))
  mn <- c(
    "Relative to SE of ref. model ELPD",
    "Relative to SE of ELPD difference"
  )
  df$method <- rep(mn, each = J + 1)
  df$method <- as.factor(df$method)
  out <- ggplot(df, aes(
    x = num_vars,
    y = value,
    color = method
  )) +
    geom_hline(yintercept = c(-1, 1), color = "gray20", lty = 2) +
    geom_line() +
    geom_point() +
    xlab("Number of terms") +
    ggtitle("ELPD relative difference") +
    theme(legend.position = "top", legend.title = element_blank())
  return(out)
}
