# MAIN DOCUMENTATION PAGE -------------------------------------------------

#' The 'lgpr2' package.
#'
#' @description Approximate longitudinal GP modeling using 'Stan'.
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Longitudinal GP Stan Bayesian
#'
#' @section Getting started:
#' See the following \code{R6} classes.
#' \itemize{
#'  \item \code{\link{TSModel}}: Main model class.
#'  \item \code{\link{TSModelFit}}: Fit class.
#'  \item \code{\link{TermList}}: Class describing model terms.
#'  \item \code{\link{FunctionDraws}}: Class to hold fitted function
#'  distributions.
#' }
#'
#' @section Data:
#' The data that you wish to analyze with 'lgpr2' should be in an \R
#' \code{data.frame} where columns correspond to measured variables and rows
#' correspond to observations. Categorical variables should be \code{factor}s
#' and continuous ones should be \code{numeric}.
#'
#' @name lgpr2-package
#' @aliases lgpr2
#' @import ggplot2 stats
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#'
#'
"_PACKAGE"


#' Run an example
#'
#' @description Fits a model to simple simulated data.
#' @export
#' @param num_bf Number of basis functions.
#' @param scale_bf Basis function domain scale.
#' @param formula The model formula.
#' @param ... Other arguments to the \code{$fit()} method of
#' \code{\link{TSModel}}.
#' @return An \code{\link{TSModelFit}} object.
example <- function(num_bf = 32, scale_bf = 1.5, formula = "y ~ gp(x)",
                    ...) {
  form <- stats::as.formula(formula)
  m <- TSModel$new(form)
  xx <- seq(1, 10, by = 0.15)
  ff <- 20 + 5 * sin(xx) + 2 * sin(5 * xx) + xx
  yy <- ff + stats::rnorm(n = length(xx), mean = 0, sd = 1)
  a <- data.frame(x = xx, y = yy)
  tc1 <- list(num_bf = num_bf, scale_bf = scale_bf)
  tc <- list(f_gp_x = tc1)
  m$fit(data = a, term_confs = tc, ...)
}
