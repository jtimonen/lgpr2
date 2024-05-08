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

#' Extend data frame and add a constant continuous variable
#'
#' @export
#' @param df original data frame
#' @param t vector of new time values
#' @param time_var name of time variable
#' @param x value of other continuous variable (will be constant)
#' @param x_var name of other continuous variable
#' @param id_var name of id variable
#' @return new data frame
extend_df2 <- function(df, t, time_var, x, x_var, id_var) {
  checkmate::assert_number(x)
  test_dat <- df
  t_test <- t
  N <- length(t_test)
  test_dat[[x_var]] <- x
  split_data <- split(test_dat, test_dat[[id_var]])
  first_rows <- lapply(split_data, function(x) x[1, ])
  test_dat <- do.call(rbind, first_rows)
  R <- nrow(test_dat)
  test_dat <- test_dat[rep(1:R, each = N), ]
  test_dat[[time_var]] <- rep(t_test, R)
  test_dat
}

#' Add input for 'FormulaSFTerm'
#'
#' @description
#' Duplicates original data columns, adding new columns with the \code{_kg} and
#' \code{_ks} suffixes if these are missing.
#'
#' @export
#' @param df the data frame
#' @param model an object of class \code{\link{TSModel}}
#' @param debug print debugging messages?
add_sff_input <- function(df, model, debug = FALSE) {
  t <- get_sff_term(model)
  vars_ks <- as.vector(t$term_list_ks$input_vars())
  vars_kg <- as.vector(t$term_list_kg$input_vars())
  vars <- c(vars_kg, vars_ks)
  for (v in vars) {
    if (v %in% colnames(df)) {
      message("column ", hl_string(v), " already exists in the data frame")
    } else {
      parts <- strsplit(v, split = "_", fixed = TRUE)[[1]]
      v_base <- paste(parts[1:(length(parts) - 1)], sep = "_")
      new_col <- df[[v_base]]
      if (is.null(new_col)) {
        stop(hl_string(v_base), " not found in the data frame")
      }
      msg <- paste0(
        "copying column ", hl_string(v_base), " to new column ",
        hl_string(v)
      )
      if (debug) {
        message(msg)
      }
      df[[v]] <- new_col
    }
  }
  df
}

#' Create a dense grid until event time for each subject (for 'JointModel')
#'
#' @export
#' @param df_tte time-to-event data frame
#' @param df_lon the longitudinal data frame
#' @param id_var id variable
#' @param time_var time variable
#' @param num_points number of grid points
#' @param even_spacing space grid points evenly?
create_jm_grid <- function(df_tte, df_lon, num_points = 30, id_var = "id",
                           time_var = "t", even_spacing = FALSE) {
  t_min <- 1e-6 # because hazard might not be defined at t = 0
  id <- df_tte[[id_var]]
  time <- df_tte[[time_var]]
  checkmate::assert_factor(id)
  checkmate::assert_numeric(time)
  ids <- as.numeric(levels(id))
  df_grid <- NULL
  t_max <- max(time)
  if (even_spacing) {
    t_seq <- seq(0, t_max, length.out = num_points)
  } else {
    t_seq <- exp(seq(0, 3, length.out = num_points)) - 1
    t_seq <- t_max / max(t_seq) * t_seq
    t_seq[t_seq <= t_min] <- t_min
  }
  for (lev in ids) {
    idx_tte <- which(id == lev)
    idx_tte <- idx_tte[length(idx_tte)]
    idx_lon <- which(df_lon[[id_var]] == lev)[1]
    row_tte <- df_tte[idx_tte, , drop = FALSE] # take last row of the subject
    row_lon <- df_lon[idx_lon, , drop = FALSE] # take first row of the subject
    row_lon[[time_var]] <- NULL # remove old time variable
    row <- cbind(row_tte, row_lon)
    df_j <- row[rep(1, num_points), ] # repeat row max_num_points times
    df_j[[time_var]] <- t_seq
    df_grid <- rbind(df_grid, df_j)
  }
  df_grid
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
