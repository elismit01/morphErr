#' Ratio Analysis Functions for morphErr
#'
#' This file contains functions for analysing ratios between measurements:
#' - mean.rcnorm.approx(): Approximates mean of ratio distribution
#' - sd.rcnorm.approx(): Approximates standard deviation of ratio distribution
#' - drcnorm.approx(): Approximation for ratio density
#' - drcnorm(): Density function for correlated normal ratios
#' - calc.mean.ratios(): Calculates mean ratios between dimensions
#'
#' @name ratio
#' @keywords internal
NULL

#' Approximate Mean of Ratio Distribution
#'
#' Calculates an approximation of the mean for a ratio of correlated normal variables.
#'
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Approximate mean of the ratio distribution
#' @keywords internal
mean.rcnorm.approx <- function(mean1, mean2, sd1, sd2, rho) {
  s <- rho * sd1/sd2
  h <- sd1 * sqrt(1 - rho^2)
  h <- smooth.indicator(mean1 - s*mean2, 0) * h
  r <- sd2/h
  b <- mean2/sd2
  a <- (mean1 - s*mean2)/h
  a/(1.01*b - 0.2713)/r + s
}

#' Approximate Standard Deviation of Ratio Distribution
#'
#' Calculates an approximation of the standard deviation for a ratio
#' of correlated normal variables.
#'
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Approximate standard deviation of the ratio distribution
#' @keywords internal
sd.rcnorm.approx <- function(mean1, mean2, sd1, sd2, rho) {
  s <- rho * sd1/sd2
  h <- sd1 * sqrt(1 - rho^2)
  h <- smooth.indicator(mean1 - s*mean2, 0) * h
  r <- sd2/h
  b <- mean2/sd2
  a <- (mean1 - s*mean2)/h
  mu <- a/(1.01*b - 0.2713)
  sqrt(((a^2 + 1)/(b^2 + 0.108*b - 3.795) - mu^2)/r^2)
}

# -------------------------------------------------------------------------------------------------------

#' Approximate Density Function for Ratio Distribution
#'
#' Calculates an approximation of the probability density function for
#' a ratio of correlated normal variables.
#'
#' @param w Vector of points at which to evaluate the density
#' @param mean1,mean2 Means of the numerator and denominator
#' @param sd1,sd2 Standard deviations of the numerator and denominator
#' @param rho Correlation between numerator and denominator
#'
#' @return Vector of density values
#' @keywords internal
drcnorm.approx <- function(w, mean1, mean2, sd1, sd2, rho) {
  c <- sd1*sd2*rho
  term1 <- dnorm((w*mean2 - mean1)/sqrt(sd1^2 - 2*w*c + w^2*sd2^2))
  term2 <- (mean2*sd1^2 - c*mean1 + (mean1*sd2^2 - c*mean2)*w)/
    ((sd1^2 - 2*w*c + w^2*sd2^2)^(3/2))
  term1*term2
}

# -------------------------------------------------------------------------------------------------------

#' Density Function for the Ratio of Correlated Normal Variables
#'
#' @description
#' Calculates the density for the ratio of two correlated normally distributed
#' random variables.
#'
#' @param w Numeric vector of values where density is to be evaluated
#' @param mean1 Mean of the numerator variable
#' @param mean2 Mean of the denominator variable
#' @param sd1 Standard deviation of the numerator variable
#' @param sd2 Standard deviation of the denominator variable
#' @param rho Correlation coefficient between numerator and denominator
#'
#' @return A numeric vector of density values corresponding to the input values
#'
#' @examples
#' w <- seq(0.5, 1.5, length.out = 10)
#' dens <- drcnorm(w, mean1 = 1, mean2 = 1, sd1 = 0.2, sd2 = 0.2, rho = 0.5)
#'
#' @export
drcnorm <- function(w, mean1, mean2, sd1, sd2, rho) {
  # Input validation
  if (!is.numeric(c(w, mean1, mean2, sd1, sd2, rho))) {
    stop("All arguments must be numeric")
  }
  if (sd1 <= 0 || sd2 <= 0) {
    stop("Standard deviations must be positive")
  }
  if (abs(rho) >= 1) {
    stop("Correlation coefficient must be between -1 and 1")
  }

  a <- (1/(1 - rho^2))^0.5*(mean1/sd1 - rho*mean2/sd2)
  b <- mean2/sd2
  tw <- (1/(1 - rho^2))^0.5*(sd2/sd1*w - rho)
  q <- (b + a*tw)/(1 + tw^2)^0.5
  phi <- dnorm(q)
  Phi <- pnorm(q) - 0.5
  f <- 1/(pi*(1 + tw^2))*exp(-0.5*(a^2 + b^2))*(1 + q/phi*Phi)
  out <- (sd2/sd1)*(1/(1 - rho^2)^0.5)*f
  out
}

# -------------------------------------------------------------------------------------------------------

#' Calculate Mean Ratios Between Dimensions
#'
#' @param fit An object of class "lme.morph" returned by fit.morph()
#' @param vcov Logical; if TRUE, returns variance-covariance matrix
#'   for ratio estimates. Default is FALSE.
#'
#' @return If vcov = FALSE, a matrix with estimates and standard errors.
#'   If vcov = TRUE, a list with components:
#'   \item{est}{Vector of ratio estimates}
#'   \item{varcov}{Variance-covariance matrix}
#'
#' @keywords internal
calc.mean.ratios <- function(fit, vcov = FALSE){
  # Get estimates w/ vcov.lme.morph
  vcov_obj <- vcov.lme.morph(fit)
  est <- vcov_obj$est
  pars.varcov <- vcov_obj$varcov

  # Extract params
  mus <- est[substr(names(est), 1, 2) == "mu"]
  sigmas <- est[substr(names(est), 1, 5) == "sigma"]
  rhos <- est[substr(names(est), 1, 3) == "rho"]
  m <- length(mus)

  # Get variances and build corrrelation matrix
  sigma.mat <- matrix(0, nrow = m, ncol = m)

  # Fill diagonal w/variances
  for(i in 1:m) {
    sigma.mat[i,i] <- sigmas[i]^2
  }

  # fill off-diagonal w/covariances
  k <- 1
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      sigma.mat[i,j] <- sigma.mat[j,i] <- rhos[k] * sigmas[i] * sigmas[j]
      k <- k + 1
    }
  }

  # Get correlation matrix
  cor.mat <- cov2cor(sigma.mat)

  # Initialise vectors for results
  n.ratios <- m * (m-1)  # Total number of ratios
  mean.ratios.est <- numeric(n.ratios)
  mean.ratios.se <- numeric(n.ratios)
  mean.ratios.names <- character(n.ratios)

  # Calculate ratios by numerator groups (each group must multiply to 1)
  ratio_groups <- list(
    # First group: dim1 as numerator (pos 1,2)
    list(
      # dim1/dim2:
      list(pos = 1, i = 1, j = 2),
      # dim1/dim3:
      list(pos = 2, i = 1, j = 3)
    ),
    # Second group: dim2 as numerator (pos3,4)
    list(
      # dim2/dim1:
      list(pos = 3, i = 2, j = 1),
      # dim2/dim3:
      list(pos = 4, i = 2, j = 3)
    ),
    # third group: dim3 as numerator (pos 5,6)
    list(
      # dim3/dim1:
      list(pos = 5, i = 3, j = 1),
      # dim3/dim2:
      list(pos = 6, i = 3, j = 2)
    )
  )

  # Processin each group to ensure ratios multiply to 1
  for(group in ratio_groups) {
    i1 <- group[[1]]$i
    j1 <- group[[1]]$j
    pos1 <- group[[1]]$pos

    ratio1 <- mean.rcnorm.approx(
      mean1 = mus[i1],
      mean2 = mus[j1],
      sd1 = sigmas[i1],
      sd2 = sigmas[j1],
      rho = cor.mat[i1,j1]
    )

    se1 <- sqrt(mus[i1]^2/mus[j1]^4 * sigma.mat[j1,j1] +
                  1/mus[j1]^2 * sigma.mat[i1,i1] -
                  2*mus[i1]/mus[j1]^3 * sigma.mat[i1,j1])

    # Calculateing second ratio to ensure product is 1
    i2 <- group[[2]]$i
    j2 <- group[[2]]$j
    pos2 <- group[[2]]$pos

    ratio2 <- 1/ratio1

    # Calculate SE for second ratio
    se2 <- sqrt(mus[i2]^2/mus[j2]^4 * sigma.mat[j2,j2] +
                  1/mus[j2]^2 * sigma.mat[i2,i2] -
                  2*mus[i2]/mus[j2]^3 * sigma.mat[i2,j2])

    # Storing results
    mean.ratios.est[pos1] <- ratio1
    mean.ratios.se[pos1] <- se1
    mean.ratios.names[pos1] <- paste0("mean(dim", i1, "/dim", j1, ")")

    mean.ratios.est[pos2] <- ratio2
    mean.ratios.se[pos2] <- se2
    mean.ratios.names[pos2] <- paste0("mean(dim", i2, "/dim", j2, ")")
  }

  # Make variance-covariance matrix for rasios
  mean.ratios.varcov <- diag(mean.ratios.se^2)

  if(vcov) {
    list(est = mean.ratios.est, varcov = mean.ratios.varcov)
  } else {
    out <- cbind(mean.ratios.est, mean.ratios.se)
    colnames(out) <- c("Estimate", "Std. Error")
    rownames(out) <- mean.ratios.names
    out
  }
}

# -------------------------------------------------------------------------------------------------------

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
