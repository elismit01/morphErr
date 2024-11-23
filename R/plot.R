#' Visualisation Functions for morphErr
#'
#' This file contains functions for visualising morphometric data and model results:
#' - plot.morph(): Plots raw morphometric data
#' - plot.lme.morph(): Creates model-based plots
#' - plot.ratio.pdf(): Plots probability density functions for ratios
#'
#' @name plot
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL

