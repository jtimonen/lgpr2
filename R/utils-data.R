#' Extend data frame
#'
#' @export
#' @param df original data frame
#' @param t vector of new time values
#' @param time_var name of time variable
#' @return new data frame
extend_df <- function(df, t, time_var) {
  u <- df_unique_factor_rows(df)
  df <- df_replicate_rows(u, length(t))
  df[[time_var]] <- rep(t, nrow(u))
  df
}


# Union of data frames
df_union <- function(x, y) {
  checkmate::assert_class(x, "data.frame")
  checkmate::assert_class(y, "data.frame")
  checkmate::assert_true(nrow(x) == nrow(y))
  common_names <- intersect(names(x), names(y))
  for (name in common_names) {
    if (!all(x[[name]] == y[[name]], na.rm = TRUE)) {
      stop(paste0("Columns '", name, "' are not identical!"))
    }
  }
  rem <- !(names(y) %in% common_names)
  rem_names <- colnames(y)[rem]
  all_names <- c(colnames(x), rem_names)
  combined_df <- cbind(x, y[, rem])
  colnames(combined_df) <- all_names
  return(combined_df)
}

# Indices of data frame columns that are factors
df_factor_cols <- function(df) {
  which(sapply(df, inherits, "factor"))
}

# Only unique rows of the factors of a data frame
df_unique_factor_rows <- function(df) {
  col_inds <- df_factor_cols(df)
  if (is.null(col_inds)) {
    stop("no factor columns found in data frame")
  }
  df <- unique(df[, col_inds, drop = FALSE])
  rownames(df) <- NULL
  df
}

# Replicate df rows
df_replicate_rows <- function(df, n) {
  return(df[rep(seq_len(nrow(df)), each = n), , drop = FALSE])
}

# Base R data frame row filter
df_filter_rows <- function(df, filter_by = NULL, kept_vals = NULL) {
  if (!is.null(filter_by)) {
    rows_keep <- which(df[[filter_by]] %in% kept_vals)
    df <- df[rows_keep, , drop = FALSE]
  }
  df
}

# Sample N data rows from each group
sample_rows_by_group <- function(data, N = 1, group_var = "arm") {
  checkmate::assert_integerish(N, lower = 1)
  group_fac <- data[[group_var]]
  checkmate::assert_factor(group_fac)
  levs <- unique(levels(group_fac))
  new_dat <- NULL
  for (lev in levs) {
    inds <- which(group_fac == lev)
    new_rows <- pick_random_rows(data[inds, , drop = F], size = N)
    new_dat <- rbind(new_dat, new_rows)
  }
  new_dat
}


# Pick random row of a data frame
pick_random_rows <- function(df, size = 1) {
  idx <- sample(nrow(df), size = size, replace = FALSE)
  df[idx, , drop = FALSE]
}

# automatic selection of prediction time points
t_pred_auto <- function(t_pred = NULL, t_data = NULL) {
  if (!is.null(t_pred)) {
    if (length(t_pred) == 1) {
      if (t_pred == "auto") {
        t_max <- min(2 * 365, max(t_data))
        message(
          "t_pred was 'auto', setting it to 10 evenly spaced points with",
          " t_max = ", t_max
        )
        t_pred <- seq(0, t_max, length.out = 10)
      } else {
        stop("invalid t_pred argument")
      }
    } else {
      checkmate::assert_numeric(t_pred, min.len = 2)
    }
  } else {
    message("t_pred was NULL, setting it based on data")
    t_range <- range(t_data)
    t_min <- min(0, t_range[1])
    t_max <- 1.1 * t_range[2]
    t_pred <- seq(t_min, t_max, length.out = 40)
  }
  t_pred
}

# Finalizer
finalize_stan_data <- function(stan_data, prior_only) {
  stan_data <- stan_data[unique(names(stan_data))] # remove duplicates
  checkmate::assert_logical(prior_only)
  stan_data$prior_only <- as.numeric(prior_only)
  stan_data
}
