# Create a new rvar with same shape as rv
rvar_set <- function(rv, value) {
  dims <- dim(posterior::as_draws_array(rv))
  posterior::rvar(x = array(value, dim = dims), dim = dim(rv))
}

# Create a new rvar with draws from standard normal, same shape as rv
rvar_std_normal <- function(rv) {
  dims <- dim(posterior::as_draws_array(rv))
  x <- stats::rnorm(n = prod(dims))
  posterior::rvar(x = array(x, dim = dims), dim = dim(rv))
}

#' Reset values of given parameters in 'rvars' list to a constant value
#'
#' @export
#' @param rvars An object of class \code{draws_rvars}.
#' @param names A character vector of parameter names.
#' @param value Value for the parameters.
#' @return An updated object of class \code{draws_rvars}.
new_draws <- function(rvars, names, value = 0) {
  checkmate::assert_character(names)
  checkmate::assert_class(rvars, "draws_rvars")
  checkmate::assert_number(value)
  for (name in names) {
    rv_orig <- rvars[[name]]
    if (is.null(rv_orig)) {
      stop("invalid parameter name '", name, "'")
    }
    message("setting ", hl_string(name), " to ", value)
    rvars[[name]] <- rvar_set(rv_orig, value)
  }
  rvars
}
