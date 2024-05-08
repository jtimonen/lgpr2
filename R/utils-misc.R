#' Name of a term in the formula syntax to its Stan code name
#'
#' @export
#' @param term A string.
term_to_code <- function(term) {
  checkmate::assert_character(term)
  form <- as.formula(paste0("y~", term))
  a <- create_termlist(form, NULL)
  a$terms[[1]]$stanname_base()
}

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
