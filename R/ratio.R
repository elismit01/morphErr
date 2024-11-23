#' Ratio Analysis Functions for morphErr
#'
#' This file contains functions for analysing ratios between measurements:
#' - calc.mean.ratios(): Calculates mean ratios between dimensions
#' - drcnorm(): Density function for correlated normal ratios
#' - drcnorm.approx(): Approximation for ratio density
#' - mean.rcnorm.approx(): Approximates mean of ratio distribution
#'
#' @name ratio
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL
