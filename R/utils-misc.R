#' Create a 'Stan' model with only given functions in it
#'
#' @param stan_file Path to the Stan file relative to the directory
#' \code{inst/functions}, without the \code{.stan} suffix.
#' @param ... Optional arguments to \code{cmdstanr::cmdstan_model()}.
stanmodel_functions <- function(stan_file, ...) {
  checkmate::assert_character(stan_file)
  sc <- read_stan_function(stan_file)
  sc <- paste("functions {", sc, "}", sep = "\n")
  sf <- cmdstanr::write_stan_file(sc)
  cmdstanr::cmdstan_model(sf, ...)
}


# Main class name
class_name <- function(x) {
  class(x)[1]
}
