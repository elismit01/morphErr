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

  # Extract parameters
  mus <- est[substr(names(est), 1, 2) == "mu"]
  m <- length(mus)

  # Get variances and build corrrelation matrix
  sigma.mat <- matrix(0, nrow = m, ncol = m)

  # Fill diagonal w/variances
  for(i in 1:m) {
    var_name <- paste0("Tau.animal.id.animal.id.var(dim", i, ")")
    sigma.mat[i,i] <- est[var_name]
  }

  # fill off-diagonal w/covariances
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      cov_name <- paste0("Tau.animal.id.animal.id.cov(dim", j, ",dim", i, ")")
      sigma.mat[i,j] <- sigma.mat[j,i] <- est[cov_name]
    }
  }

  # Standard deviations and correlations
  sigmas <- sqrt(diag(sigma.mat))
  cor.mat <- cov2cor(sigma.mat)

  # Calculate ratios + prepare for ses
  n.ratios <- 2*choose(m, 2)
  mean.ratios.est <- numeric(n.ratios)
  mean.ratios.names <- character(n.ratios)
  mean.ratios.se <- numeric(n.ratios)

  # Ratios and se using delta methode
  k <- 1
  mean.ratios.varcov <- matrix(0, n.ratios, n.ratios)

  for(i in 1:m) {
    for(j in 1:m) {
      if(i != j) {
        # Mean ratio
        mean.ratios.est[k] <- mean.rcnorm.approx(
          mean1 = mus[i],
          mean2 = mus[j],
          sd1 = sigmas[i],
          sd2 = sigmas[j],
          rho = cor.mat[i,j]
        )

        # Approx se using delta method
        mean.ratios.se[k] <- sqrt(mus[i]^2/mus[j]^4 * sigma.mat[j,j] +
                                    1/mus[j]^2 * sigma.mat[i,i] -
                                    2*mus[i]/mus[j]^3 * sigma.mat[i,j])

        mean.ratios.names[k] <- paste0("mean(dim", i, "/dim", j, ")")
        k <- k + 1
      }
    }
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

