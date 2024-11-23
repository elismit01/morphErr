#' Prediction Functions for morphErr
#'
#' This file contains functions for making predictions from fitted models:
#' - predict.lme.morph(): Predicts measurements from fitted models
#' - calc.betas(): Calculates coefficients for predictions
#'
#' @name predict
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL
