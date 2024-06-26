# Hilbert space basis functions
hs_bf <- function(x, b, L) {
  k <- (pi * b) / (2 * L)
  1 / sqrt(L) * sin(k * (x + L))
}

# Wiggliness penalty
hs_wiggliness_diag <- function(m, L) {
  a <- (pi * m) / (2 * L)
  a^4
}

# Matrix S for HS GP (with c=1)
hs_penalty_matrix <- function(B, L) {
  S_diag <- hs_wiggliness_diag(seq_len(B), L)
  n <- length(S_diag)
  diag(S_diag, nrow = n, ncol = n)
}

# Get L
get_domain_size <- function(object) {
  L <- object$xt$dom_size
  if (is.null(L)) {
    stop("Domain size L not given!")
  }
  checkmate::assert_number(L, lower = 0)
  L
}


#' @importFrom mgcv smooth.construct
#' @export
mgcv::smooth.construct

#' @importFrom mgcv Predict.matrix
#' @export
mgcv::Predict.matrix

#' Hilbert space GP smooth
#' @export
smooth.construct.hs.smooth.spec <- function(object, data, knots) {
  ## object$p.order = null space dimension

  if (object$bs.dim < 0) object$bs.dim <- 10 ## default
  x <- data[[object$term]] ## the data
  B <- object$bs.dim
  L <- get_domain_size(object)
  if (B <= 1) {
    stop("Must have at least 2 basis functions")
  }

  X <- matrix(0, length(x), B)
  for (b in seq_len(B)) {
    X[, b] <- hs_bf(x, b, L)
  }
  object$X <- X # the finished model matrix
  if (!object$fixed) {
    # create the penalty matrix
    S <- hs_penalty_matrix(B, L)
    object$S[[1]] <- S
  }
  object$rank <- B # penalty rank
  # object$null.space.dim <- B - 1 # dim. of unpenalized space
  # object$df <- ncol(object$X) # maximum DoF (if unconstrained)

  class(object) <- "hs.smooth" # Give object a class
  object
}

#' prediction method function for the `hs' smooth class
#' @export
Predict.matrix.hs.smooth <- function(object, data) {
  x <- data[[object$term]]
  B <- object$bs.dim
  X <- matrix(0, length(x), B)
  L <- get_domain_size(object)
  for (b in seq_len(B)) {
    X[, b] <- hs_bf(x, b, L)
  }
  X # return the prediction matrix
}

# HS smooth by
create_hs_smooth <- function(x, z = NULL, L = 1.5, B = 32) {
  bs_str <- "hs"
  long <- paste0(", bs='", bs_str, "', k=", B, ", xt=list(dom_size=", L, ")")
  by <- paste0(", by = ", z)
  if (is.null(z)) {
    by <- ""
  }
  paste0("s(", x, long, by, ")")
}
