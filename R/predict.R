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

  # Estimates if nt provided
  if (is.null(est)) {
    vcov_obj <- vcov.lme.morph(fit)
    est <- vcov_obj$est
  }

  # Get parameters using new format
  mus <- est[substr(names(est), 1, 2) == "mu"]
  sigmas <- est[substr(names(est), 1, 5) == "sigma"]
  rhos <- est[substr(names(est), 1, 3) == "rho"]

  # Constructing cov matrix
  m <- length(mus)
  sigma.mat <- matrix(0, nrow = m, ncol = m)

  # Filll diagonal w/variances
  for(i in 1:m) {
    sigma.mat[i,i] <- sigmas[i]^2
  }

  # Fill off-diagonal w/covariances
  k <- 1
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      sigma.mat[i,j] <- sigma.mat[j,i] <- rhos[k] * sigmas[i] * sigmas[j]
      k <- k + 1
    }
  }

  if (type == "lm") {
    # lm interpretation
    sub.sigma.mat.y <- sigma.mat[y.dim, x.dim, drop = FALSE]
    sub.sigma.mat.x <- sigma.mat[x.dim, x.dim, drop = FALSE]
    beta.dims <- c(sub.sigma.mat.y %*% solve(sub.sigma.mat.x))
    beta0 <- mus[y.dim] - sum(beta.dims * mus[x.dim])
    betas.est <- c(beta0, beta.dims)
  } else {
    # PCA/RMA interpretation (to fix error it now uses sigmas directly)
    betas.est <- c(mus[y.dim] - sigmas[y.dim]/sigmas[x.dim] * mus[x.dim],
                   sigmas[y.dim]/sigmas[x.dim])
  }

  # Handle ses if requested
  if (stders) {
    # Delta method for se
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

# -------------------------------------------------------------------------------------------------------

#' Predict measurements using true measurements
#'
#' @description
#' Makes predictions for one dimension based on true measurements of other dimensions.
#'
#' @param object A fitted model object from fit.morph()
#' @param y.dim Integer specifying which dimension to predict
#' @param newdata A data frame of other dimensions to use for prediction. Column names
#'        must be of the form "dimX" where X is the dimension number
#' @param type Either "lm" or "pca" for prediction type
#' @param ... Additional arguments passed to methods
#'
#' @return A matrix with Estimate and Std.Error columns
#' @export
predict.lme.morph <- function(object, y.dim, newdata = NULL, type = c("lm", "pca"), ...) {
  # Input validation
  if (!inherits(object, "lme.morph")) {
    stop("Invalid model object. Input must be a fitted model of class 'lme.morph'")
  }

  # Valid types
  valid_types <- c("lm", "pca")
  if (!type[1] %in% valid_types) {
    stop(
      "Invalid type argument. See ?predict.lme.morph for possible selections."
    )
  }
  type <- match.arg(type)

  # If no measurements, return pop means with se
  if (is.null(newdata)) {
    vcov.obj <- object$vcov
    est <- vcov.obj$est
    mus <- est[substr(names(est), 1, 2) == "mu"]
    pred.est <- mus[y.dim]
    pred.se <- sqrt(diag(vcov.obj$varcov))[paste0("mu", y.dim)]

    result <- matrix(c(pred.est, pred.se), nrow = 1)
    colnames(result) <- c("Estimate", "Std. Error")
    return(result)
  }

  # Validate newdata
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame with columns named 'dimX' where X is the dimension number")
  }

  # Get predictor dimes
  x.dim <- as.numeric(sapply(strsplit(colnames(newdata), "dim"), function(x) x[2]))

  # Get coeficients and vcov matrix
  betas <- calc.betas(fit = object, y.dim = y.dim, x.dim = x.dim, type = type, vcov = FALSE)
  vcov.obj <- calc.betas(fit = object, y.dim = y.dim, x.dim = x.dim, type = type, vcov = TRUE)
  vcov.mat <- vcov.obj$varcov

  # Create pred matrix and calculate prediction
  X <- as.matrix(cbind(1, newdata))
  pred.est <- sum(betas[,"Estimate"] * c(1, unlist(newdata)))

  # Calculate se w/delta method
  X.mat <- matrix(c(1, unlist(newdata)), nrow = 1)
  pred.se <- sqrt(X.mat %*% vcov.mat %*% t(X.mat))

  result <- matrix(c(pred.est, pred.se), nrow = 1)
  colnames(result) <- c("Estimate", "Std. Error")
  return(result)
}

# -------------------------------------------------------------------------------------------------------

#' Predict measurements from observed measurements with error
#'
#' @description
#' Makes predictions by combining prior information with observed measurements that contain error.
#'
#' @param object A fitted model object from fit.morph()
#' @param true A vector of true dimension measurements. Set elements to NA for dimensions to predict
#' @param obs A matrix of observed measurements with error, where each row represents measurements
#'        from one photograph. NA values allowed for missing measurements.
#'
#' @return A numeric vector of predictions for all dimensions
#' @export
predict.from.obs <- function(object, true, obs = NULL) {
  # Input validation
  if (!inherits(object, "lme.morph")) {
    stop("Invalid model object. Input must be a fitted model of class 'lme.morph'")
  }

  # Get params
  vcov.obj <- object$vcov
  est <- vcov.obj$est
  m <- length(levels(object$data$dim))
  pars <- organise.pars(est, n.animals = 1, n.photos = 1, m = m, block.only = TRUE)

  # Validate tru vector length
  if (length(true) != m) {
    stop(
      "Length of 'true' must match number of dimensions in model"
    )
  }

  # Handle observations
  if (!is.null(obs)) {
    if (is.vector(obs)) {
      obs <- matrix(obs, nrow = 1)
    }
    if (ncol(obs) != m) {
      stop(
        "Number of columns in 'obs' must match number of dimensions in model"
      )
    }
  } else {
    obs <- matrix(nrow = 0, ncol = m)
  }

  # Determine which dims to predict
  dims.to.predict <- which(is.na(true))
  if (length(dims.to.predict) == 0) {
    stop("No dimensions to predict (no NA values in 'true' vector)")
  }

  # Optimisation function
  djoint <- function(t.to.predict) {
    t <- true
    t[dims.to.predict] <- t.to.predict
    out <- 0

    # Add observation likelihood if ther are obs
    if (nrow(obs) > 0) {
      for (j in 1:nrow(obs)) {
        obs.dims <- !is.na(obs[j, ])
        out <- out + mvtnorm::dmvnorm(obs[j, obs.dims],
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
  preds <- nlminb(pars$mus[dims.to.predict], djoint)

  # Prepare output
  # Start with true vals:
  out <- true
  # Fill in preds:
  out[dims.to.predict] <- preds$par

  return(out)
}
