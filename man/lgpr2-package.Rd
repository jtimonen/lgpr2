% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa.R
\docType{package}
\name{lgpr2-package}
\alias{lgpr2-package}
\alias{lgpr2}
\title{The 'lgpr2' package.}
\description{
Approximate longitudinal GP modeling using 'Stan'.
}
\section{Getting started}{

See the following \code{R6} classes.
\itemize{
 \item \code{\link{LonModel}}: Main model class.
 \item \code{\link{LonModelFit}}: Fit class.
 \item \code{\link{TermList}}: Class describing model terms.
 \item \code{\link{FunctionDraws}}: Class to hold fitted function
 distributions.
}
}

\section{Data}{

The data that you wish to analyze with 'lgpr2' should be in an \R
\code{data.frame} where columns correspond to measured variables and rows
correspond to observations. Categorical variables should be \code{factor}s
and continuous ones should be \code{numeric}.
}

\author{
Juho Timonen (first.last at iki.fi)
}
\keyword{Bayesian}
\keyword{GP}
\keyword{Longitudinal}
\keyword{Stan}
