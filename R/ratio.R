#' Ratio Analysis Functions for morphErr
#'
#' This file contains functions for analysing ratios between measurements:
#' - calc.mean.ratios(): Calculates mean ratios between dimensions
#'
#' @name ratio
#' @keywords internal
NULL

#' Calculate Conditional Ratios Between Dimensions
#'
#' @description
#' Calculates a ratio between dimensions, conditional on the denominator.
#' Used internally by plotting functions.
#'
#' @param fit A fitted model object of class "lme.morph"
#' @param y.dim Integer specifying numerator dimension
#' @param x.dim Integer specifying denominator dimension
#' @param newdata.x.dim Numeric vector of values for denominator dimension
#' @param type Character string, either "lm" or "pca", specifying calculation method
#'
#' @return A matrix with columns:
#'   \item{Estimate}{Conditional ratio estimates}
#'   \item{Std. Error}{Standard errors for the estimates}
#'
#' @details
#' Calculates ratios between dimensions conditional on the denominator dimension value.
#' Can use either linear model or PCA interpretation for the calculations.
#'
#' @keywords internal
calc.conditional.ratio <- function(fit, y.dim, x.dim, newdata.x.dim, type = "lm") {
  # Input validation
  if (!inherits(fit, "lme.morph")) {
    stop("'fit' must be an object of class 'lme.morph'")
  }

  # Check dims
  data <- getData(fit)
  m <- length(unique(data$dim))
  if (x.dim > m || y.dim > m) {
    stop("Dimension index exceeds number of dimensions in data")
  }
  if (x.dim == y.dim) {
    stop("Numerator and denominator dimensions must be different")
  }

  # Extract param estimates
  est <- summary(fit)[, 1]
  mus <- est[substr(names(est), 1, 2) == "mu"]
  sigmas <- est[substr(names(est), 1, 5) == "sigma"]
  rhos <- est[substr(names(est), 1, 3) == "rho"]
  psis <- est[substr(names(est), 1, 3) == "psi"]
  phis <- est[substr(names(est), 1, 3) == "phi"]
  par.list <- list(mus = mus, sigmas = sigmas, rhos = rhos, psis = psis, phis = phis)
  par.vec <- unlist(par.list)

  # Initialise gradient matric
  conditional.ratios.gr <- matrix(0, nrow = length(newdata.x.dim), ncol = length(par.vec))

  if (type == "lm") {
    # lm interpretation
    rho.name <- paste0("rho", sort(c(x.dim, y.dim))[1], ",", sort(c(x.dim, y.dim))[2])
    rho <- rhos[rho.name]

    # Calculate estimates
    conditional.ratios.est <- mus[y.dim]/newdata.x.dim +
      rho*sigmas[y.dim]/sigmas[x.dim]*(1 - mus[x.dim]/newdata.x.dim)

    # Calculate gradints
    conditional.ratios.gr[, names(est) == paste0("mu", y.dim)] <- 1/newdata.x.dim
    conditional.ratios.gr[, names(est) == paste0("mu", x.dim)] <-
      -rho*sigmas[y.dim]/(sigmas[x.dim]*newdata.x.dim)
    conditional.ratios.gr[, names(est) == paste0("sigma", y.dim)] <-
      rho/sigmas[x.dim] - rho*mus[x.dim]/(sigmas[x.dim]*newdata.x.dim)
    conditional.ratios.gr[, names(est) == paste0("sigma", x.dim)] <-
      -rho*sigmas[y.dim]/(sigmas[x.dim]^2) +
      rho*sigmas[y.dim]*mus[x.dim]/(sigmas[x.dim]^2*newdata.x.dim)
    conditional.ratios.gr[, names(est) == rho.name] <-
      sigmas[y.dim]/sigmas[x.dim] -
      sigmas[y.dim]*mus[x.dim]/(sigmas[x.dim]*newdata.x.dim)

  } else if (type == "pca") {
    # PCA interpretation
    conditional.ratios.est <- mus[y.dim]/newdata.x.dim +
      sigmas[y.dim]/sigmas[x.dim]*(1 - mus[x.dim]/newdata.x.dim)

    # Calculate gradients
    conditional.ratios.gr[, names(est) == paste0("mu", y.dim)] <- 1/newdata.x.dim
    conditional.ratios.gr[, names(est) == paste0("mu", x.dim)] <-
      -sigmas[y.dim]/(sigmas[x.dim]*newdata.x.dim)
    conditional.ratios.gr[, names(est) == paste0("sigma", y.dim)] <-
      1/sigmas[x.dim] - mus[x.dim]/(sigmas[x.dim]*newdata.x.dim)
    conditional.ratios.gr[, names(est) == paste0("sigma", x.dim)] <-
      -sigmas[y.dim]/(sigmas[x.dim]^2) +
      sigmas[y.dim]*mus[x.dim]/(sigmas[x.dim]^2*newdata.x.dim)
  }

  # Calculate ses using delta method
  pars.varcov <- fit$vcov$varcov
  conditional.ratios.varcov <- conditional.ratios.gr %*% pars.varcov %*% t(conditional.ratios.gr)
  conditional.ratios.se <- sqrt(diag(conditional.ratios.varcov))

  # Prep output
  out <- cbind(conditional.ratios.est, conditional.ratios.se)
  colnames(out) <- c("Estimate", "Std. Error")
  out
}
