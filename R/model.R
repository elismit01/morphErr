#' Model Fitting Functions for morphErr
#'
#' This file contains the main model fitting functions for morphometric analysis:
#' - fit.morph(): Main function for fitting morphometric models
#' - organise.pars(): Organises parameter vectors into lists
#' - construct.varcov(): Constructs variance-covariance matrices
#'
#' @name model
#' @keywords internal
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl
#' @importFrom lmeInfo Fisher_info extract_varcomp
NULL

#' Organise Parameter Vector into Components
#'
#' Takes a vector of parameters and organises it into separate objects
#' (e.g., a vector of mu parameters, variance-covariance matrices
#' sigma and xi, etc).
#'
#' @param pars A numeric vector of parameters
#' @param n.animals Integer, number of animals
#' @param n.photos Integer vector, number of photos per animal
#' @param m Integer, number of dimensions measured
#' @param block.only Logical, if TRUE returns only the block matrix
#'
#' @return A list containing:
#'   \item{mus}{Vector of mean parameters}
#'   \item{sigma}{Variance-covariance matrix for animal effects}
#'   \item{xi}{Variance-covariance matrix for measurement error}
#'
#' @keywords internal
#'
organise.pars <- function(pars, n.animals, n.photos, m, block.only = FALSE) {
  n.pars <- length(pars)
  n.mus <- m
  n.sigmas <- m
  n.rhos <- sum(1:(m - 1))
  n.psis <- m
  n.phis <- sum(1:(m - 1))
  if (n.pars != sum(n.mus + n.sigmas + n.rhos + n.psis + n.phis)){
    stop("Incorrect number of elements in pars.")
  }
  mus <- numeric(n.mus)
  k <- 1
  for (i in 1:n.mus){
    mus[i] <- pars[k]
    k <- k + 1
  }
  sigmas <- numeric(n.sigmas)
  for (i in 1:n.sigmas){
    sigmas[i] <- pars[k]
    k <- k + 1
  }
  rhos <- numeric(n.rhos)
  for (i in 1:n.rhos){
    rhos[i] <- pars[k]
    k <- k + 1
  }
  psis <- numeric(n.psis)
  for (i in 1:n.psis){
    psis[i] <- pars[k]
    k <- k + 1
  }
  phis <- numeric(n.phis)
  for (i in 1:n.phis){
    phis[i] <- pars[k]
    k <- k + 1
  }
  sigma.mat <- construct.varcov(sigmas, rhos, n.animals, m, block.only)
  xi.mat <- construct.varcov(psis, phis, sum(n.photos), m, block.only)
  list(mus = mus, sigma = sigma.mat, xi = xi.mat)
}

#' Construct Variance-Covariance Matrix
#'
#' Takes standard deviations and correlation parameters, and organises
#' them into a block-diagonal variance-covariance matrix.
#'
#' @param sds Numeric vector of standard deviations for each dimension
#' @param cors Numeric vector of correlation parameters between dimensions
#' @param n.blocks Integer, number of blocks in the matrix
#' @param m Integer, number of dimensions measured
#' @param block.only Logical, if TRUE returns only one block rather than the full matrix
#'
#' @return A variance-covariance matrix
#'
#' @keywords internal
#'
construct.varcov <- function(sds, cors, n.blocks, m, block.only){
  block <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m){
    block[i, i] <- sds[i]^2
  }
  k <- 1
  for (i in 1:(m - 1)){
    for (j in (i + 1):m){
      block[i, j] <- block[j, i] <- cors[k]*sds[i]*sds[j]
      k <- k + 1
    }
  }
  if (block.only){
    out <- block
  } else {
    out <- matrix(0, nrow = m*n.blocks, ncol = m*n.blocks)
    for (i in 1:n.blocks){
      out[((i - 1)*m + 1):((i - 1)*m + m),
          ((i - 1)*m + 1):((i - 1)*m + m)] <- block
    }
  }
  out
}
