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

# Base R data frame row filter
df_filter_rows <- function(df, filter_by = NULL, kept_vals = NULL) {
  if (!is.null(filter_by)) {
    rows_keep <- which(df[[filter_by]] %in% kept_vals)
    df <- df[rows_keep, , drop = FALSE]
  }
  df
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

# Finalizer
finalize_stan_data <- function(stan_data, prior_only) {
  stan_data <- stan_data[unique(names(stan_data))] # remove duplicates
  checkmate::assert_logical(prior_only)
  stan_data$prior_only <- as.numeric(prior_only)
  stan_data
}
