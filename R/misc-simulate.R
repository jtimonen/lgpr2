#' Data simulation for model reduction experiments
#'
#' @export
#' @param n number of data points
#' @param rho correlation
#' @param rel_true true relevances
#' @param sigma noise magnitude
simulate <- function(n, rho, rel_true, sigma) {
  D <- length(rel_true) / 2
  checkmate::assert_integerish(D)
  x <- simulate_x(n, rho, D)
  # D <- length(rel_true)
  # x <- matrix(runif(n = n * D), n, D)
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  f <- simulate_f(x)
  f <- decide_relevant(f, rel_true)
  f_tot <- rowSums(f)
  y <- f_tot + sigma * rnorm(n = length(f_tot))
  list(
    x = x,
    f = f,
    f_tot = f_tot,
    y = y
  )
}



# Add unrelated variables
create_correlated_x <- function(x, p) {
  s <- stats::sd(x)
  rnorm(n = length(x), mean = x, sd = p * s)
}

# Flip categorical
flip_z_of_id <- function(z, id, id_to_flip) {
  levs <- as.numeric(levels(z))
  ids <- as.numeric(levels(id))
  rows_edit <- which(as.numeric(id) %in% id_to_flip)
  cat_to_flip <- unique(as.numeric(z)[rows_edit])
  hmm <- setdiff(levs, cat_to_flip)
  new_cat <- sample(as.character(hmm), size = 1)
  z_new <- as.numeric(z)
  z_new[rows_edit] <- as.numeric(new_cat)
  as.factor(z_new)
}

# Add unrelated variables
create_correlated_z <- function(id, z, n_edit) {
  levs <- as.numeric(levels(z))
  ids <- as.numeric(levels(id))
  ids_remaining <- ids
  for (n in seq_len(n_edit)) {
    id_to_flip <- sample(ids_remaining, size = 1)
    z <- flip_z_of_id(z, id, id_to_flip)
    ids_remaining <- setdiff(ids_remaining, id_to_flip)
  }
  z
}

# Add unrelated variables
create_new_x <- function(bn, dat, p) {
  for (j in 1:length(p)) {
    x_base <- dat[[bn]]
    xn <- paste0(bn, "_u", j)
    dat[[xn]] <- create_correlated_x(x_base, p = p[j])
  }
  dat
}

# Add unrelated variables
create_new_z <- function(bn, dat, n) {
  for (j in 1:length(n)) {
    z_base <- dat[[bn]]
    zn <- paste0(bn, "_u", j)
    dat[[zn]] <- create_correlated_z(dat$id, z_base, n[j])
  }
  dat
}

# Add unrelated variables
add_unrelated <- function(dat, n_x = 4, n_z = 4, n_flip_z = 4) {
  p_x <- seq(0.75, 1.0, length.out = n_x)
  n_z_edit <- rep(n_flip_z, n_z)
  dat <- create_new_x("x", dat, p_x)
  # dat <- create_new_x("x2", dat, p_x)
  dat <- create_new_z("z", dat, n_z_edit)
  # dat <- create_new_z("z2", dat, n_z_edit)
  cnd <- colnames(dat)
  list(
    dat = dat,
    xn = cnd[which(grepl(cnd, pattern = "x"))],
    zn = cnd[which(grepl(cnd, pattern = "z"))],
    px = p_x,
    nze = n_z_edit
  )
}

# Covariates
simulate_x <- function(n, rho, D = 4) {
  Sigma <- rho * matrix(1, D, D) + (1 - rho) * diag(nrow = D, ncol = D)
  mu <- rep(0, nrow(Sigma))
  x1 <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  x2 <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  cbind(x1, x2)
}

# Function components
simulate_f <- function(x) {
  f <- 0 * x
  for (j in 1:ncol(f)) {
    a <- 0.2 + 1.2 * runif(1)
    b <- pi * runif(1)
    f_j <- sin(a * x[, j] + b)
    f_j <- f_j - mean(f_j)
    f_j <- f_j / stats::sd(f_j)
    f[, j] <- f_j
  }
  f
}

# Multiply columns by relevance
decide_relevant <- function(f, rel) {
  N <- dim(f)[1]
  D <- dim(f)[2]
  f * matrix(rel, nrow = N, ncol = D, byrow = TRUE)
}


# Create data split
create_split <- function(n, p_test) {
  checkmate::assert_integerish(n, lower = 1)
  checkmate::assert_number(p_test, lower = 0, upper = 1)
  p_train <- 1 - p_test
  n_train <- round(n * p_train)
  idx_train <- sort(sample(n, n_train))
  list(
    train = idx_train,
    test = sort(setdiff(seq_len(n), idx_train))
  )
}

# Take a subset of data
subset_data <- function(sim, inds) {
  list(
    x = sim$x[inds, ],
    f = sim$f[inds, ],
    f_tot = sim$f_tot[inds],
    y = sim$y[inds]
  )
}

# Create data frame for other methods
create_data_frame <- function(ds, test = FALSE) {
  df <- cbind(ds$dat$x, ds$dat$y)
  df <- as.data.frame(df)
  idx <- ds$split$train
  if (test) {
    idx <- ds$split$test
  }
  df <- df[idx, ]
  D <- dim(df)[2] - 1
  colnames(df) <- c(paste("x", seq_len(D), sep = ""), "y")
  df
}
