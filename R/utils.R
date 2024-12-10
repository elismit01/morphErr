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

#' Constructor Function for Mean Ratio Calculations
#'
#' @param dim1,dim2 Integers specifying which dimensions to use
#'
#' @return A function that calculates ratios
#' @keywords internal
mean.ratio.constructor <- function(dim1, dim2) {
  function(pars) {
    # Get parameters
    mus <- pars$mus
    sigmas <- pars$sigmas
    rhos <- pars$rhos

    # Plain matrix operations not RTMB
    m <- length(mus)
    sigma.mat <- matrix(0, nrow = m, ncol = m)

    # Fill diagonal
    for (i in 1:m) {
      sigma.mat[i,i] <- sigmas[i]^2
    }

    # Fill corelations
    k <- 1
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
        sigma.mat[i,j] <- sigma.mat[j,i] <- rhos[k] * sigmas[i] * sigmas[j]
        k <- k + 1
      }
    }

    # Get correlation matrix and calculate mean ratio
    cor.mat <- cov2cor(sigma.mat)
    mean.rcnorm.approx(mus[dim1], mus[dim2],
                       sqrt(sigma.mat[dim1,dim1]),
                       sqrt(sigma.mat[dim2,dim2]),
                       cor.mat[dim1,dim2])
  }
}


#' Extract Variance-Covariance Matrix from Morphometric Model
#'
#' S3 method for extracting both parameter estimates and their
#' variance-covariance matrix from a fitted morphometric model.
#'
#' @param object A fitted model of class "lme.morph"
#'
#' @return A list containing:
#'   \item{est}{Named vector of parameter estimates}
#'   \item{varcov}{Variance-covariance matrix}
#'
#' @keywords internal
vcov.lme.morph <- function(object) {
  # Convert lme.morph objectinto standard lme
  object.lme <- object
  class(object.lme) <- "lme"

  # Get original parameter estimates
  orig.est <- c(object$coefficients$fixed,
                extract_varcomp(object, vector = TRUE))
  names(orig.est)[1:length(object$coefficients$fixed)] <-
    paste0("mu", 1:length(object$coefficients$fixed))

  # Getting original variance-covariance matrix
  orig.varcov <- as.matrix(bdiag(vcov(object.lme),
                                 solve(Fisher_info(object))))

  # Get Components and create list
  list(est = orig.est,
       varcov = orig.varcov)
}
