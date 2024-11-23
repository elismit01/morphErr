#' Simulation Functions for morphErr
#'
#' This file contains functions for simulating morphometric data:
#' - sim.measurements(): Simulates measurement data from a survey
#' - sim.morph(): Simulates multiple datasets and fits models
#' - extract.sim.morph(): Extracts results from simulation studies
#'
#' @name simulate
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL

