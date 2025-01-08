#' Prediction Functions for morphErr
#'
#' This file contains functions for making predictions from fitted models:
#' - calc.betas(): Calculates coefficients for predictions
#' - predict.lme.morph(): Predicts measurements from fitted models
#'
#' @name predict
#' @keywords internal
NULL

#' Calculate Beta Coefficients for Predictions
#'
#' @description
#' Calculates coefficients for predicting one dimension from others,
#' using either linear model or PCA interpretation.
#'
#' @param fit A fitted model object of class "lme.morph"
#' @param est Optional parameter estimates (if NULL, extracted from fit)
#' @param stders Logical; whether to compute standard errors
#' @param y.dim Integer specifying which dimension to predict
#' @param x.dim Integer vector specifying which dimensions to use as predictors
#' @param type Character string, either "lm" or "pca"
#' @param vcov Logical; if TRUE, returns variance-covariance matrix
#'
#' @return If vcov = FALSE, a matrix with estimates and standard errors.
#'   If vcov = TRUE, a list with components:
#'   \item{est}{Vector of coefficient estimates}
#'   \item{varcov}{Variance-covariance matrix}
#'
#' @keywords internal
calc.betas <- function(fit, est = NULL, stders = TRUE, y.dim, x.dim,
                       type = "lm", vcov = FALSE) {
  # Check type for PCA case
  if (type == "pca" && length(x.dim) > 1) {
    stop("Type 'pca' is only available for single predictor relationships")
  }

  # Get estimates if not provided
  if (is.null(est)) {
    vcov_obj <- vcov.lme.morph(fit)
    est <- vcov_obj$est
  }

  # Get parameters
  mus <- est[substr(names(est), 1, 2) == "mu"]
  sigmas <- sqrt(est[grep("var", names(est))])

  # Extract corelations from covariances
  var_names <- names(est)[grep("var", names(est))]
  cov_names <- names(est)[grep("cov", names(est))]
  rhos <- numeric(length(cov_names))
  for(i in seq_along(cov_names)) {
    dims <- as.numeric(regmatches(cov_names[i],
                                  regexpr("\\d+", cov_names[i])))
    rhos[i] <- est[cov_names[i]] /
      (sqrt(est[var_names[dims[1]]]) * sqrt(est[var_names[dims[2]]]))
  }

  # Construct cov matrix
  m <- length(mus)
  sigma.mat <- matrix(0, nrow = m, ncol = m)
  for(i in 1:m) {
    sigma.mat[i,i] <- sigmas[i]^2
  }
  k <- 1
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      sigma.mat[i,j] <- sigma.mat[j,i] <- rhos[k] * sigmas[i] * sigmas[j]
      k <- k + 1
    }
  }

  if (type == "lm") {
    # Linear model interpretation
    sub.sigma.mat.y <- sigma.mat[y.dim, x.dim, drop = FALSE]
    sub.sigma.mat.x <- sigma.mat[x.dim, x.dim, drop = FALSE]
    beta.dims <- c(sub.sigma.mat.y %*% solve(sub.sigma.mat.x))
    beta0 <- mus[y.dim] - sum(beta.dims * mus[x.dim])
    betas.est <- c(beta0, beta.dims)
  } else {
    # PCA/RMA interpretation
    betas.est <- c(mus[y.dim] - sigmas[y.dim]/sigmas[x.dim] * mus[x.dim],
                   sigmas[y.dim]/sigmas[x.dim])
  }

  # Handle se's if requested
  if (stders) {
    # Use delta method for se (simplified approximation)
    if (type == "lm") {
      se0 <- sqrt(sigma.mat[y.dim, y.dim])
      se.dims <- sqrt(diag(solve(sub.sigma.mat.x)) * sigma.mat[y.dim, y.dim])
      betas.se <- c(se0, se.dims)
    } else {
      se0 <- sqrt(sigma.mat[y.dim, y.dim] +
                    (sigmas[y.dim]/sigmas[x.dim])^2 * sigma.mat[x.dim, x.dim])
      se1 <- sigmas[y.dim]/sigmas[x.dim] *
        sqrt(sigma.mat[y.dim, y.dim]/sigmas[y.dim]^2 +
               sigma.mat[x.dim, x.dim]/sigmas[x.dim]^2)
      betas.se <- c(se0, se1)
    }

    if (vcov) {
      # Create approx vcov matrix
      vcov.mat <- diag(betas.se^2)
      out <- list(est = betas.est, varcov = vcov.mat)
    } else {
      out <- cbind(betas.est, betas.se)
      colnames(out) <- c("Estimate", "Std. Error")
      rownames(out) <- c("beta0", paste0("beta", x.dim))
    }
  } else {
    out <- betas.est
  }

  out
}

#' Predict measurements for morphometric data
#'
#' @description
#' Makes predictions based on either true measurements or observed measurements with error.
#' Users must provide either true_measurements or observed_measurements, not both.
#'
#' @param object A fitted model object from fit.morph()
#' @param true_measurements A data frame of true measurements to predict from, with column names
#'        of the form "dimX" where X is the dimension number, or NULL
#' @param observed_measurements A matrix of observed measurements with error, where each row
#'        represents measurements from one photograph, or NULL. NA values allowed for missing
#'        measurements.
#' @param y.dim Which dimension to predict
#' @param type Either "lm" or "pca" for prediction type (only used with true measurements)
#' @param ... Additional arguments passed to methods
#'
#' @return A matrix with Estimate and Std.Error columns
predict.lme.morph <- function(object,
                              true_measurements = NULL,
                              observed_measurements = NULL,
                              y.dim,
                              type = c("lm", "pca"),
                              ...) {
  # Input validation
  if (!inherits(object, "lme.morph")) {
    stop("Object must be of class 'lme.morph'")
  }

  if (!is.null(true_measurements) && !is.null(observed_measurements)) {
    stop("Cannot provide both true and observed measurements. Please provide only one type.")
  }

  type <- match.arg(type)

  # Handle true measurements path
  if (!is.null(true_measurements)) {
    if (!is.data.frame(true_measurements)) {
      stop("true_measurements must be a data frame")
    }

    x.dim <- as.numeric(sapply(strsplit(colnames(true_measurements), "dim"), function(x) x[2]))

    # Get coefficients (not vcov matrix)
    betas <- calc.betas(fit = object, y.dim = y.dim, x.dim = x.dim, type = type, vcov = FALSE)

    # Create prediction matrix
    X <- as.matrix(cbind(1, true_measurements))

    # Calculate prediction and standard error
    pred.est <- sum(betas[,"Estimate"] * c(1, unlist(true_measurements)))
    pred.se <- NA  # We'll implement SE calculation later if needed

    result <- matrix(c(pred.est, pred.se), nrow = 1)
    colnames(result) <- c("Estimate", "Std. Error")
    return(result)
  }

  # Handle observed measurements path
  if (!is.null(observed_measurements)) {
    if (is.vector(observed_measurements)) {
      observed_measurements <- matrix(observed_measurements, nrow = 1)
    }
    if (is.data.frame(observed_measurements)) {
      observed_measurements <- as.matrix(observed_measurements)
    }

    # Get parameter vals
    vcov.obj <- object$vcov
    est <- vcov.obj$est
    m <- length(levels(object$data$dim))
    pars <- organise.pars(est, n.animals = 1, n.photos = 1, m = m, block.only = TRUE)

    # Create prediction vector
    true <- numeric(m)
    # All dimensions to be predicted:
    true[] <- NA
    # Fill non-predicted dimensions with means:
    true[-y.dim] <- pars$mus[-y.dim]

    # Setup optimisation function
    djoint <- function(t.to.predict) {
      t <- true
      t[y.dim] <- t.to.predict
      out <- 0

      # Add observation likelihood
      if (nrow(observed_measurements) > 0) {
        for (j in 1:nrow(observed_measurements)) {
          obs.dims <- !is.na(observed_measurements[j, ])
          out <- out + mvtnorm::dmvnorm(observed_measurements[j, obs.dims],
                                        mean = t[obs.dims],
                                        sigma = pars$xi[obs.dims, obs.dims, drop = FALSE],
                                        log = TRUE)
        }
      }

      # Add prior
      out <- out + mvtnorm::dmvnorm(t, mean = pars$mus, sigma = pars$sigma, log = TRUE)
      -out
    }

    # Optimise
    preds <- nlminb(pars$mus[y.dim], djoint)
    pred.est <- preds$par

    # Return with NA for se
    result <- matrix(c(pred.est, NA), nrow = 1)
    colnames(result) <- c("Estimate", "Std. Error")
    return(result)
  }

  # If no measurements provided, return pop means with SE
  vcov.obj <- object$vcov
  est <- vcov.obj$est
  mus <- est[substr(names(est), 1, 2) == "mu"]
  pred.est <- mus[y.dim]
  pred.se <- sqrt(diag(vcov.obj$varcov))[paste0("mu", y.dim)]

  result <- matrix(c(pred.est, pred.se), nrow = 1)
  colnames(result) <- c("Estimate", "Std. Error")
  return(result)
}

