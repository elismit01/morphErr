% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{calc.betas}
\alias{calc.betas}
\title{Calculate Beta Coefficients for Predictions}
\usage{
calc.betas(
  fit,
  est = NULL,
  stders = TRUE,
  y.dim,
  x.dim,
  type = "lm",
  vcov = FALSE
)
}
\arguments{
\item{fit}{A fitted model object of class "lme.morph"}

\item{est}{Optional parameter estimates (if NULL, extracted from fit)}

\item{stders}{Logical; whether to compute standard errors}

\item{y.dim}{Integer specifying which dimension to predict}

\item{x.dim}{Integer vector specifying which dimensions to use as predictors}

\item{type}{Character string, either "lm" or "pca"}

\item{vcov}{Logical; if TRUE, returns variance-covariance matrix}
}
\value{
If vcov = FALSE, a matrix with estimates and standard errors.
If vcov = TRUE, a list with components:
\item{est}{Vector of coefficient estimates}
\item{varcov}{Variance-covariance matrix}
}
\description{
Calculates coefficients for predicting one dimension from others,
using either linear model or PCA interpretation.
}
\keyword{internal}
