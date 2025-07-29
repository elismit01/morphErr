# Utility Functions for morphErr
#
# This file contains helper functions and utilities:
# - print.summary.lme.morph(): Prints model summaries
# - summary.lme.morph(): Creates model summaries
# - smooth.indicator(): Helper function for smooth approximations


# -------------------------------------------------------------------------------------------------------

#' Extract Variance-Covariance Matrix from Morphometric Model
#'
#' An S3 method for extracting both parameter estimates and their
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
    data <- getData(object)
    ## Number of dimensions.
    m <- length(unique(data$dim))
    ## Extracting the estimtaes, but the parameterisation is
    ## different from what we want: (1) we have elements of the
    ## variance-covariance matrix of the random effects, instead of
    ## standard deviations and correlations; and (2) we have an
    ## estimate of the measurement error standard deviation for the
    ## first dimension with multipliers for the other dimensions,
    ## instead of having separate standard deviations for each
    ## dimension, and (3) if we haven't fitted an intercept, then we
    ## need to adjust the betas to provide the mus.
    orig.est <- c(object$coefficients$fixed, extract_varcomp(object, vector = TRUE))
    names(orig.est)[1:m] <- paste0("beta", 1:m)
    orig.names <- names(orig.est)
    ## Creating a vector to hold all parameter estimates.
    n.est <- length(orig.est)
    est <- numeric(n.est)
    dim.pairs <- apply(combn(m, 2), 2, function(x) paste(x, collapse = ","))
    names(est) <- c(paste0("mu", 1:m), paste0("sigma", 1:m), paste0("rho", dim.pairs),
                    paste0("psi", 1:m), paste0("phi", dim.pairs))
    ## Creating a list with functions to transform these estimates to our
    ## parameterisation.
    g <- vector("list", length = n.est)
    k <- 1
    ## Filling in est and g for mu parameters.
    for (i in 1:m){
        orig.name.beta <- paste0("beta", i)
        which.orig.name.beta <- which(orig.names == orig.name.beta)
        if (object$intercept & i != 1){
            which.orig.name.intercept <- which(orig.names == "beta1")
            est[k] <- orig.est["beta1"] + orig.est[orig.name.beta]
            g[[k]] <- reformulate(paste0("x", which.orig.name.intercept, "+ x", which.orig.name.beta))
        } else {
            est[k] <- orig.est[orig.name.beta]
            g[[k]] <- reformulate(paste0("x", which.orig.name.beta))
        }
        k <- k + 1
    }
    ## Filling in est and g for sigma parameters.
    for (i in 1:m){
        orig.name <- paste0("Tau.animal.id.animal.id.var(dim", i, ")")
        which.orig.name <- which(orig.names == orig.name)
        est[k] <- sqrt(orig.est[orig.name])
        g[[k]] <- reformulate(paste0("sqrt(x", which.orig.name, ")"))
        k <- k + 1
    }
    ## Filling in est and g for rho parameters.
    for (i in 1:(m - 1)){
        for (j in (i + 1):m){
            orig.name.sd1 <- paste0("Tau.animal.id.animal.id.var(dim", i, ")")
            orig.name.sd2 <- paste0("Tau.animal.id.animal.id.var(dim", j, ")")
            orig.name.cov <- paste0("Tau.animal.id.animal.id.cov(dim", j, ",dim", i, ")")
            which.orig.name.sd1 <- which(orig.names == orig.name.sd1)
            which.orig.name.sd2 <- which(orig.names == orig.name.sd2)
            which.orig.name.cov <- which(orig.names == orig.name.cov)
            est[k] <- orig.est[orig.name.cov]/(sqrt(orig.est[orig.name.sd1])
                *sqrt(orig.est[orig.name.sd2]))
            g[[k]] <- reformulate(paste0("x", which.orig.name.cov, "/(sqrt(x",
                                         which.orig.name.sd1, ")*sqrt(x", which.orig.name.sd2,
                                         "))"))
            k <- k + 1
        }
    }
    ## Filling in est and g for psi parameters.
    orig.name.ss <- "sigma_sq"
    which.orig.name.ss <- which(orig.names == "sigma_sq")
    for (i in 1:m){
        if (i == 1){
            est[k] <- sqrt(orig.est[orig.name.ss])
            g[[k]] <- reformulate(paste0("sqrt(x", which.orig.name.ss, ")"))
        } else {
            orig.name.var <- paste0("var_params", i - 1)
            which.orig.name.var <- which(orig.names == orig.name.var)
            est[k] <- sqrt(orig.est[orig.name.ss])*orig.est[orig.name.var]
            g[[k]] <- reformulate(paste0("sqrt(x", which.orig.name.ss, ")*x", which.orig.name.var))
        }
        k <- k + 1
    }
    ## Filling in est and g for phi parameters.
    l <- 1
    for (i in 1:(m - 1)){
        for (j in (i + 1):m){
            orig.name.cor <- paste0("cor_params", l)
            which.orig.name.cor <- which(orig.names == orig.name.cor)
            est[k] <- orig.est[orig.name.cor]
            g[[k]] <- reformulate(paste0("x", which.orig.name.cor))
            l <- l + 1
            k <- k + 1
        }
    }
    object.lme <- object
    class(object.lme) <- "lme"
    ## Original variance-covariance matrix. It's block diagonal because
    ## fixed-effect estimators are independent of variance component
    ## estimators.
    orig.varcov <- as.matrix(bdiag(vcov(object.lme), solve(Fisher_info(object))))
    ## Obtaining variance-covariance matrix under our parameterisation.
    varcov <- deltamethod(g, orig.est, orig.varcov, ses = FALSE)
    rownames(varcov) <- colnames(varcov) <- names(est)
    list(est = est, varcov = varcov)
}

# -------------------------------------------------------------------------------------------------------

#' @export
print.summary.lme.morph <- function(x, ...) {
  class(x) <- "matrix"
  x <- as.data.frame(x)
  if (any(names(x) == "P-value")) {
    x["P-value"] <- format.pval(x["P-value"])
  }
  print(x)
}

# -------------------------------------------------------------------------------------------------------

#' Summarise Morphometric Model Fits
#'
#' @description An S3 method that creates summaries of fitted
#'     morphometric models. This function can provide different types
#'     of summaries including parameter estimates, beta coefficients,
#'     and tests for isometry.
#'
#' @details The following summaries are available with this method,
#'     controlled by the `type` argument:
#'
#' \describe{
#' 
#' \item{`type = "pars"`}{Estimates and standard errors for the model
#'                        parameters, including the vectors \eqn{\mu}
#'                        for the means of the true dimension sizes,
#'                        \eqn{\sigma} for the standard deviations for
#'                        true dimension sizes, \eqn{\rho} for the
#'                        correlations between true dimension sizes,
#'                        \eqn{\psi} for the standard deviations in
#'                        measurement errors for the dimensions, and
#'                        \eqn{\phi} for the correlations between
#'                        measurement errors for all pairs of
#'                        dimensions.}
#' 
#' \item{`type = "betas"` or `type = "betas-lm"`}{Estimates and
#'                                                standard errors for
#'                                                coefficients of a
#'                                                linear combination
#'                                                that provides the
#'                                                expected value of
#'                                                one dimension
#'                                                (specified by
#'                                                `y.dim`) using any
#'                                                subset of the
#'                                                remaining dimensions
#'                                                (specified in the
#'                                                vector `x.dim`).}
#'
#' \item{`type = "betas-pca"`}{An estimated intercept and slope for
#'                             the reduced major axis (or principal
#'                             component axis) summarising the
#'                             relationship between the true values of
#'                             dimensions specified by `y.dim` and
#'                             `x.dim`.}
#'
#' \item{`type = "isometric-pca"`}{Results for tests of the null
#'                                 hypotheses of isometric
#'                                 relationships between all pairs of
#'                                 dimensions. The p-values are
#'                                 obtained using normal
#'                                 approximations for the test
#'                                 statistics.}
#'
#' \item{`type = "isometric-pca-boot"`}{The same as
#'                                      `"isometric-pca"`, except a
#'                                      bootstrap is used to obtain
#'                                      p-values. The normal
#'                                      approximation and bootstrap
#'                                      typically provide very similar
#'                                      results.}
#' 
#' }
#'
#' @param object An object of class `lme.morph`, returned by
#'     [`fit.morph()`].
#' @param ... Other parameters (for S3 generic compatibility).
#' @param type A character string specifying type of summary. Either
#'     `"pars"` (the default), `"betas"` `"betas-lm"`, `"betas-pca"`,
#'     "`isometric-pca`", or `"isometric-pca-boot"`. See 'Details' for
#'     further information.
#' @param y.dim An integer specifying the response dimension for when
#'     `type` is `"betas"`, `"betas-lm"`, or `"betas-pca"`.
#' @param x.dim An integer vector specifying the explanatory
#'     dimensions for when `type` is `"betas"`, `"betas-lm"`, or
#'     `"betas-pca"`.
#' @param B Number of bootstrap iterations for when `type` is
#'     `"isometric-pca-boot"`.
#'
#' @return A matrix or data frame containing the requested summary.
#' @examples
#' ## Fitting model to manta ray data.
#' fit <- fit.morph(manta)
#' ## Parameter estimates and standard errors.
#' summary(fit)
#' ## Estimated coefficients for the linear combination that provides
#' ## the expected value of the first dimension from only the second.
#' summary(fit, type = "betas-lm", y.dim = 1, x.dim = 2)
#' ## Estimated coefficients for the linear combination that provides
#' ## the expected value of the second dimension from both the first
#' ## and third.
#' summary(fit, type = "betas-lm", y.dim = 2, x.dim = c(1, 3))
#' ## Estimated intercept and slope for the reduced major axis (or
#' ## pricipal component axis) summarising the relationship between
#' ## the first and second dimensions.
#' summary(fit, type = "betas-pca", y.dim = 1, x.dim = 2)
#' ## Tests for isometry between all pairs of dimensions.
#' summary(fit, type = "isometric-pca")
#' @export
summary.lme.morph <- function(object, ..., type = "pars", y.dim, x.dim, B = 1000){
  # Valid types
  valid_types <- c(
    "pars",   # Param estimates
    "betas", "betas-lm",  # lm coefficients
    "betas-pca",   # PCA interpretaton
    "isometric-pca",  # Isometry test results
    "isometric-pca-boot"  # Bootstrapped isometry test results
  )

  # Check if type is valid
  if (!type %in% valid_types) {
    stop(
      "Invalid type argument. See ?summary.lme.morph for possible selections."
    )
  }

  ## Indicator for whether the model involved a log transformation.
  log.transform <- object$log.transform
    
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
          out.names[k] <- paste0("dim", i, " vs dim", j)
          if (type == "isometric-pca-boot") {        
            betas.boot <- apply(boots, 1, function(x)
                calc.betas(est = x, stders = FALSE, y.dim = j,
                           x.dim = i, type = "pca"))
            if (log.transform) {
              betas.boot <- betas.boot[2, ]
              out[k, ] <- calc.betas(fit = object, est = est, y.dim = j, x.dim = i, type = "pca")[2, ]        
              p.val.onesided <- mean(betas.boot >= 1)
              p.val.onesided <- min(c(p.val.onesided, 1 - p.val.onesided))
              p.val[k] <- 2*p.val.onesided
            } else {
              betas.boot <- betas.boot[1, ]
              out[k, ] <- calc.betas(fit = object, est = est, y.dim = j, x.dim = i, type = "pca")[1, ]  
              p.val.onesided <- mean(betas.boot >= 0)
              p.val.onesided <- min(c(p.val.onesided, 1 - p.val.onesided))
              p.val[k] <- 2*p.val.onesided
            }
          } else if (type == "isometric-pca") {
            if (log.transform){
              out[k, ] <- calc.betas(fit = object, est = est, y.dim = j, x.dim = i, type = "pca")[2, ]
              z <- (out[k, 1] - 1)/out[k, 2]
              p.val[k] <- 2*pnorm(-abs(z))
            } else {
              ## Using the distribution of alpha0.
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










