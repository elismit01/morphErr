# Model Fitting Functions for morphErr
#
# This file contains the main model fitting functions for morphometric analysis:
# - fit.morph(): Main function for fitting morphometric models
# - organise.pars(): Organises parameter vectors into lists
# - construct.varcov(): Constructs variance-covariance matrices


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

# -------------------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------------------

#' Fit Morphometric Model
#'
#' Fits the model described by \href{../doc/model.pdf}{Stevenson,
#' Smit, and Setyawan (in submission)} to morphometric data.
#'
#' @details
#' 
#' This is a special case of a linear mixed-effects model, and is
#' fitted via a call to [`nlme::lme()`]. Some arguments of
#' `fit.morph()` are directly passed to [`nlme::lme()`].
#' 
#' @param method A character string indicating the objective function
#'     used to fit the model. Either `"ML"` for maximum likelihood or
#'     `"REML"` for restricted maximum likelihood.
#' @param control A list of control values for the estimation
#'     algorithm to replace the default values returned by the
#'     function [`nlme::lmeControl()`].
#'
#' @return An object of classes `lme.morph` (specific to this package)
#'     and `lme` (inherited from [`nlme::lme()`]). The object is the
#'     list returned by [`nlme::lme()`] with a couple of additional
#'     components.
#'
#' Rather than inspecting the object directly, the best way to extract
#' understandable output from the object is using the S3 methods
#' [`summary.lme.morph()`], [`plot.lme.morph()`], and
#' [`predict.lme.morph()`] via the generic functions `summary()`,
#' `plot()`, and `predict()`. Other S3 methods for the `lme` class are
#' availble via `nlme`, such as [`nlme::coef.lme()`].
#'
#' @examples
#' \dontrun{
#' ## Fitting model to manta ray data.
#' fit <- fit.morph(manta)
#' ## Using maximum likelihood instead of REML.
#' fit <- fit.morph(manta, method = "ML")
#' }
#' @inheritSection plotmorph The `data` argument
#' @inheritParams plotmorph
#' @export
fit.morph <- function(data, method = "REML",
                      control = list(maxIter = 100000,
                                     msMaxIter = 100000)){
  # Ensure inputs are factors.
  for (i in c("animal.id", "photo.id", "dim")){
    if (!is.factor(data[, i])){
      data[, i] <- factor(data[, i])
    }
  }

  # Reordering the data because (believe it or not) the
  # parameterisation of the variance model in nlme uses whichever
  # dimension appears first in the data as the baseline.
  data <- data[order(data$dim), ]
    
  # Set up groups
  gdata <- groupedData(measurement ~ 1 | animal.id / photo.id, data = data)

  # Set up fixed effects formula
  fixed.arg <- measurement ~ 0 + dim

  # Fit model
  fit <- lme(fixed = fixed.arg,
             random = ~ 0 + dim | animal.id,
             correlation = corSymm(form = ~ 1 | animal.id / photo.id),
             weights = varIdent(form = ~ 1 | dim),
             data = gdata,
             control = control,
             method = method)

  # Add additional attributes
  fit$intercept <- FALSE
  class(fit) <- c("lme.morph", class(fit))
  fit$vcov <- vcov(fit)

  return(fit)
}
