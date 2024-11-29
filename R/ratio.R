#' Ratio Analysis Functions for morphErr
#'
#' This file contains functions for analysing ratios between measurements:
#' - mean.rcnorm.approx(): Approximates mean of ratio distribution
#' - sd.rcnorm.approx(): Approximates standard deviation of ratio distribution
#' - drcnorm.approx(): Approximation for ratio density
#' - drcnorm(): Density function for correlated normal ratios
#' - calc.mean.ratios(): Calculates mean ratios between dimensions
#'
#' @name ratio
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL

#' Approximate Mean of Ratio Distribution
#'
#' Calculates an approximation of the mean for a ratio of correlated normal variables.
#'
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Approximate mean of the ratio distribution
#' @keywords internal
mean.rcnorm.approx <- function(mean1, mean2, sd1, sd2, rho) {
  s <- rho * sd1/sd2
  h <- sd1 * sqrt(1 - rho^2)
  h <- smooth.indicator(mean1 - s*mean2, 0) * h
  r <- sd2/h
  b <- mean2/sd2
  a <- (mean1 - s*mean2)/h
  a/(1.01*b - 0.2713)/r + s
}

#' Approximate Standard Deviation of Ratio Distribution
#'
#' Calculates an approximation of the standard deviation for a ratio
#' of correlated normal variables.
#'
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Approximate standard deviation of the ratio distribution
#' @keywords internal
sd.rcnorm.approx <- function(mean1, mean2, sd1, sd2, rho) {
  s <- rho * sd1/sd2
  h <- sd1 * sqrt(1 - rho^2)
  h <- smooth.indicator(mean1 - s*mean2, 0) * h
  r <- sd2/h
  b <- mean2/sd2
  a <- (mean1 - s*mean2)/h
  mu <- a/(1.01*b - 0.2713)
  sqrt(((a^2 + 1)/(b^2 + 0.108*b - 3.795) - mu^2)/r^2)
}

#' Approximate Density Function for Ratio Distribution
#'
#' Calculates an approximation of the probability density function for
#' a ratio of correlated normal variables.
#'
#' @param w Vector of points at which to evaluate the density
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Vector of density values
#' @keywords internal
drcnorm.approx <- function(w, mean1, mean2, sd1, sd2, rho) {
  c <- sd1*sd2*rho
  term1 <- dnorm((w*mean2 - mean1)/sqrt(sd1^2 - 2*w*c + w^2*sd2^2))
  term2 <- (mean2*sd1^2 - c*mean1 + (mean1*sd2^2 - c*mean2)*w)/
    ((sd1^2 - 2*w*c + w^2*sd2^2)^(3/2))
  term1*term2
}
