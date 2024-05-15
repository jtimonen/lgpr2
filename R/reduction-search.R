# Base class
ForwardSearch <- R6Class(
  classname = "ForwardSearch",
  public = list(
    num_comps = NULL,
    num_steps = 1,
    verbosity = 0,
    score_name = "unknown",
    mlpd = list(),
    path = NULL, # possible predefined path

    # Init
    initialize = function(J, path = NULL, verbosity = 1, num_steps = 1) {
      self$num_comps <- J
      self$verbosity <- verbosity
      self$path <- path
      self$num_steps <- num_steps
    },

    # Score a model
    score = function(model, ...) {
      cat("* Inheriting class should override the score method!\n")
      return(runif(1))
    },

    # Compute score for each candidate
    step = function(model, candidates, ...) {
      scores <- rep(0, length(candidates))
      jj <- 0
      for (ca in candidates) {
        jj <- jj + 1
        c_model <- c(model, ca)
        scores[jj] <- self$score(c_model, ...)
      }

      scores
    },

    # Get candidate models at current step
    get_candidates = function(cur_model, J) {
      j <- length(cur_model) + 1
      if (!is.null(self$path)) {
        message("Using pre-defined search path")
        cands <- self$path[j] # use pre-defined path
      } else {
        message("Possibly many alternative models at this step")
        cands <- setdiff(seq_len(J), cur_model)
      }
    },

    # Perform search
    run = function(...) {
      J <- self$num_comps
      model <- NULL
      score <- NULL
      while (length(model) < J) {
        cands <- self$get_candidates(model, J)
        cand_scores <- self$step(model, cands, ...)
        i_best <- which(cand_scores == max(cand_scores))
        score <- c(score, cand_scores[i_best])
        model <- c(model, cands[i_best])
        if (self$verbosity > 0) {
          print(model)
        }
      }
      list(path = model, score = score)
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
    score = function(model, refmod, eval_mode, kl0, elpd_loo_ref) {
      res <- refmod$project(model, eval_mode = eval_mode)
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
    step = function(model, candidates, refmod) {
      if (length(candidates) == 1) {
        return(1)
      }
      df <- NULL
      elpd_loo_ref <- refmod$get_loo_estimates()
      for (ca in candidates) {
        c_model <- c(model, ca)
        res <- self$score(c_model, refmod, FALSE, NULL, elpd_loo_ref)
        df <- rbind(df, res)
      }
      df <- data.frame(df)
      cand_scores <- df$score
      which(cand_scores == max(cand_scores))
    },

    # Perform search
    run = function(refmod, ...) {
      J <- self$num_comps
      model <- NULL
      score <- NULL
      elpd <- NULL
      elpd_loo_ref <- refmod$get_loo_estimates()
      history <- self$score(model, refmod, TRUE, NULL, elpd_loo_ref)
      kl0 <- history$kl
      # thresh <- getOption("pp.threshold", default = 0.95)

      # Loop
      j <- 0
      while (length(model) < J) {
        j <- j + 1
        msg <- paste0("Step ", j, "/", J, ".")
        message(msg)

        # Choose best candidate
        cands <- self$get_candidates(model, J)
        i_best <- self$step(model, cands, refmod)

        # Score the best candidate
        new_model <- c(model, cands[i_best])
        new_row <- self$score(new_model, refmod, TRUE, kl0, elpd_loo_ref)
        print(new_row)

        # Update history
        history <- rbind(history, new_row)
        model <- c(model, cands[i_best])
        edr <- new_row$elpd_loo_rel_diff
        cat("elpd_loo_rel_diff = ", edr, "\n", sep = "")
        if (j >= self$num_steps) {
          break
        }
      }

      # Return
      list(path = model, history = history)
    }
  )
)


# Projection predictive forward search
pp_forward_search <- function(refmod, path, num_steps) {
  J <- length(refmod$get_model()$term_list$terms)
  a <- ProjectionForwardSearch$new(J, path, num_steps = num_steps)
  a$run(refmod)
}

# Plot result of above
plot_pp_pexp <- function(res, thresh = 0.95) {
  J <- length(res$path)
  num_vars <- c(0, 1:J)
  p_exp <- res$history$p_explained
  df <- data.frame(num_vars, p_exp)
  out <- ggplot(df, aes(x = num_vars, y = p_exp)) +
    geom_hline(yintercept = thresh, color = "firebrick3", lty = 2) +
    ylab("1 - KL / KL0")
  out <- out + geom_line() + geom_point() +
    xlab("Number of variables") + ggtitle("Forward search")
  return(out)
}

# Plot result of above
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
    xlab("Number of variables") +
    ggtitle("Forward search", subtitle = st)
  return(out)
}
