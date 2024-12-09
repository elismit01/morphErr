#' Utility Functions for morphErr
#'
#' This file contains helper functions and utilities:
#' - print.summary.lme.morph(): Prints model summaries
#' - summary.lme.morph(): Creates model summaries
#' - smooth.indicator(): Helper function for smooth approximations
#'
#' @name utils
#' @keywords internal
NULL

#' Smooth Indicator Function
#'
#' Creates a smoothed version of an indicator function.
#'
#' @param x Numeric value to evaluate
#' @param t Threshold value
#'
#' @return A smoothed indicator value between 0 and 1
#' @keywords internal
smooth.indicator <- function(x, t){
  2*pnorm(x - t, sd = 1e-10) - 1
}
