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

#' Fit Morphometric Model
#'
#' Fits a linear mixed-effects model to morphometric data accounting for measurement error.
#'
#' @param data A data frame containing:
#'   \itemize{
#'     \item animal.id: factor indicating which animal the measurement is from
#'     \item photo.id: factor indicating which photo the measurement is from
#'     \item dim: factor indicating which dimension is measured
#'     \item measurement: the observed morphometric measurement
#'   }
#' @param method Character, either "ML" for maximum likelihood or "REML" for
#'   restricted maximum likelihood
#' @param intercept Logical, whether to include an intercept in fixed effects
#'
#' @return An object of class "lme.morph" containing:
#'   \itemize{
#'     \item coefficients: Model coefficients
#'     \item vcov: Variance-covariance matrix
#'     \item residuals: Model residuals
#'     \item fitted: Fitted values
#'   }
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' data <- sim.measurements(n.animals = 10, n.photos = rep(5, 10), m = 3,
#'                         c(290, 130, 75, 45, 25, 15,
#'                           0.85, 0.90, 0.95,
#'                           2.00, 1.50, 1.00,
#'                           0.40, 0.50, 0.60))
#'
#' # Fit the model
#' fit <- fit.morph(data, method = "REML")
#' }
#'
#' @export
fit.morph <- function(data, method = "REML", intercept = FALSE){
  # Ensure input sare factors
  for (i in c("animal.id", "photo.id", "dim")){
    if (!is.factor(data[, i])){
      data[, i] <- factor(data[, i])
    }
  }

  # Set up groups
  gdata <- nlme::groupedData(measurement ~ 1 | animal.id / photo.id, data = data)

  # Set up fixed effects formula
  if (intercept){
    fixed.arg <- measurement ~ dim
  } else {
    fixed.arg <- measurement ~ 0 + dim
  }

  # Fit model
  fit <- nlme::lme(fixed = fixed.arg,
                   random = ~ 0 + dim | animal.id,
                   correlation = nlme::corSymm(form = ~ 1 | animal.id / photo.id),
                   weights = nlme::varIdent(form = ~ 1 | dim),
                   data = gdata,
                   control = nlme::lmeControl(maxIter = 10000, msMaxIter = 10000),
                   method = method)

  # Add additional attributes
  fit$intercept <- intercept
  class(fit) <- c("lme.morph", class(fit))
  fit$vcov <- vcov(fit)

  return(fit)
}
