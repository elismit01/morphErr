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


#' Print Method for Morphometric Model Summary
#'
#' @description
#' Prints the summary of a morphometric model in a clean format.
#'
#' @param x Object of class "summary.lme.morph"
#' @param ... Additional arguments passed to print methods
#'
#' @return No return value, called for side effect of printing
#' @export
print.summary.lme.morph <- function(x, ...) {
  class(x) <- "matrix"
  x <- as.data.frame(x)
  if (any(names(x) == "P-value")) {
    x["P-value"] <- format.pval(x["P-value"])
  }
  print(x)
}


#' Summary Method for Morphometric Models
#'
#' @description
#' Creates summaries of fitted morphometric models. Can provide different types of
#' summaries including parameter estimates, beta coefficients, and isometry tests.
#'
#' @param object A fitted model of class "lme.morph"
#' @param ... Additional arguments passed to methods
#' @param type Character string specifying type of summary:
#'   - "pars": parameter estimates and standard errors
#'   - "betas"/"betas-lm": coefficients for linear model interpretation
#'   - "betas-pca": coefficients for PCA interpretation
#'   - "ratios": mean ratios between dimensions
#'   - "isometric-pca": isometry test results
#'   - "isometric-pca-boot": bootstrapped isometry test results
#' @param y.dim Integer specifying which dimension to predict (for beta coefficients)
#' @param x.dim Integer vector specifying predictor dimensions (for beta coefficients)
#' @param B Number of bootstrap iterations for isometric-pca-boot
#' @param boot.invert Logical; whether to invert bootstrap distribution
#'
#' @return A matrix or data frame containing the requested summary statistics
#' @export
summary.lme.morph <- function(object, ..., type = "pars", y.dim, x.dim, B = 1000, boot.invert = FALSE) {
  vcov.obj <- object$vcov
  ## Get estimates
  est <- vcov.obj$est
  if (type == "pars") {
    ## Extracting SEs
    ses <- sqrt(diag(vcov.obj$varcov))
    ##  Estimates and ses for model parameters
    out <- cbind(est, ses)
    colnames(out) <- c("Estimate", "Std. Error")
  } else if (type == "betas" | type == "betas-lm") {
    ## Beta parameters for predicting y.dim from x.dims
    out <- calc.betas(fit = object, est = est, y.dim = y.dim, x.dim = x.dim)
  } else if (type == "betas-pca") {
    ## Beta parameters for the first PC between y.dim and x.dim
    out <- calc.betas(fit = object, est = est, y.dim = y.dim, x.dim = x.dim, type = "pca")
  } else if (type == "ratios") {
    ## Estimates and SES for mean ratios between dimensions
    out <- calc.mean.ratios(object)
  } else if (type == "isometric-pca" | type == "isometric-pca-boot") {
    data <- getData(object)
    ## num of dimensions
    m <- length(unique(data$dim))
    ## Tests for isometric growth between all dimensions
    n.comparisons <- 2*choose(m, 2)
    out <- matrix(0, nrow = n.comparisons, ncol = 2)
    out.names <- character(n.comparisons)
    if (type == "isometric-pca-boot") {
      boots <- rmvnorm(B, est, vcov.obj$varcov)
    }
    p.val <- numeric(n.comparisons)
    k <- 1
    for (i in 1:m) {
      for (j in 1:m) {
        if (i != j) {
          ## Extracting beta0.
          out[k, ] <- calc.betas(fit = object, est = est, y.dim = j, x.dim = i, type = "pca")[1, ]
          out.names[k] <- paste0("dim", i, " vs dim", j)
          if (type == "isometric-pca-boot") {
            betas.boot <- apply(boots, 1, function(x) calc.betas(est = x, stders = FALSE, y.dim = j, x.dim = i, type = "pca"))[1, ]
            ## The ol' flipperooni
            if (boot.invert) {
              betas.boot <- 2*out[k, 1] - betas.boot
            }
            p.val.onesided <- mean(betas.boot >= 0)
            p.val.onesided <- min(c(p.val.onesided, 1 - p.val.onesided))
            p.val[k] <- 2*p.val.onesided
          } else if (type == "isometric-pca") {
            ## Using the distribution of alpha0
            mu.mean <- vcov.obj$est[paste0("mu", c(j, i))]
            mu.varcov <- vcov.obj$varcov[paste0("mu", c(j, i)), paste0("mu", c(j, i))]
            sigma.mean <- vcov.obj$est[paste0("sigma", c(j, i))]
            sigma.varcov <- vcov.obj$varcov[paste0("sigma", c(j, i)), paste0("sigma", c(j, i))]
            ## Calculating this product thing, which
            ## has a true value of zero under
            ## isometric growth.
            out[k, 1] <- mu.mean[1]*sigma.mean[2] - mu.mean[2]*sigma.mean[1]
            ## Getting relevant parameter names and
            ## positions in the estimates vector.
            par.names <- names(vcov.obj$est)
            par.name.mu1 <- paste0("mu", j)
            which.par.name.mu1 <- which(par.names == par.name.mu1)
            par.name.mu2 <- paste0("mu", i)
            which.par.name.mu2 <- which(par.names == par.name.mu2)
            par.name.sigma1 <- paste0("sigma", j)
            which.par.name.sigma1 <- which(par.names == par.name.sigma1)
            par.name.sigma2 <- paste0("sigma", i)
            which.par.name.sigma2 <- which(par.names == par.name.sigma2)
            ## Setting up the expression for the deltamethod() function.
            g <- reformulate(paste0("x", which.par.name.mu1, " * x", which.par.name.sigma2,
                                    " - x", which.par.name.mu2, " * x", which.par.name.sigma1))
            ## Doing the delta method to get a standard error.
            out[k, 2] <- deltamethod(g, vcov.obj$est, vcov.obj$varcov)
            p.val.lower <- pnorm(out[k, 1], sd = out[k, 2])
            p.val[k] <- 2*min(c(p.val.lower, 1 - p.val.lower))
          }
          k <- k + 1
        }
      }
    }
    out <- cbind(out, p.val)
    rownames(out) <- out.names
    colnames(out) <- c("Estimate", "Std. Error", "P-value")
  }
  class(out) <- "summary.lme.morph"
  out
}










