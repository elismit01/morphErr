# Prediction Functions for morphErr
#
# This file contains functions for making predictions from fitted models:
# - predict.lme.morph(): Predicts measurements from fitted models using true measurements
# - predict.obs(): Predicts measurements from fitted models using observed measurements


# -------------------------------------------------------------------------------------------------------

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
                                        ## Check type for PCA case
    if (type == "pca" && length(x.dim) > 1) {
        stop("Type 'pca' is only available for single predictor relationships")
    }
    
    ## Estimates if nt provided
    if (is.null(est)) {
        vcov_obj <- vcov.lme.morph(fit)
        est <- vcov_obj$est
    }
    
    ## Get parameters using new format
    mus <- est[substr(names(est), 1, 2) == "mu"]
    sigmas <- est[substr(names(est), 1, 5) == "sigma"]
    rhos <- est[substr(names(est), 1, 3) == "rho"]
    psis <- est[substr(names(est), 1, 3) == "psi"]
    phis <- est[substr(names(est), 1, 3) == "phi"]
    par.list <- list(mus = mus, sigmas = sigmas, rhos = rhos, psis = psis, phis = phis)
    par.vec <- unlist(par.list)
    
    ## Constructing cov matrix
    m <- length(mus)
    sigma.mat <- matrix(0, nrow = m, ncol = m)
    
    ## Fill diagonal w/variances
    for(i in 1:m) {
        sigma.mat[i,i] <- sigmas[i]^2
    }
    
    ## Fill off-diagonal w/covariances
    k <- 1
    for(i in 1:(m-1)) {
        for(j in (i+1):m) {
            sigma.mat[i,j] <- sigma.mat[j,i] <- rhos[k] * sigmas[i] * sigmas[j]
            k <- k + 1
        }
    }

    if (type == "lm") {
        ## lm interpretation
        sub.sigma.mat.y <- sigma.mat[y.dim, x.dim, drop = FALSE]
        sub.sigma.mat.x <- sigma.mat[x.dim, x.dim, drop = FALSE]
        beta.dims <- c(sub.sigma.mat.y %*% solve(sub.sigma.mat.x))
        beta0 <- mus[y.dim] - sum(beta.dims * mus[x.dim])
        betas.est <- c(beta0, beta.dims)
    } else {
                                        ## PCA/RMA interpretation (to fix error it now uses sigmas directly)
        betas.est <- c(mus[y.dim] - sigmas[y.dim]/sigmas[x.dim] * mus[x.dim],
                       sigmas[y.dim]/sigmas[x.dim])
    }

  ## Handle ses if requested
    if (stders) {
        ## Delta method for se
        if (type == "lm") {
            n.x <- length(x.dim)
            ## Inverse of the submatrix.
            sigmaqq.inv <- solve(sub.sigma.mat.x)
            ## Matrix to store the partial derivatives.
            betas.gr <- matrix(0, nrow = n.x + 1, ncol = length(par.vec))
            ## Partial derivatives for the intercept, starting with the mus.
            mu.deriv <- numeric(length(mus))
            mu.deriv[y.dim] <- 1
            mu.deriv[x.dim] <- -sub.sigma.mat.y %*% sigmaqq.inv
            ## Now the sigmas.
            sigma.deriv <- numeric(length(sigmas))
            dsigmapq <- matrix(0, nrow = 1, ncol = n.x)
            for (i in 1:n.x){
                dsigmapq[1, i] <- get.corpar(rhos, y.dim, x.dim[i])*sigmas[x.dim[i]]
            }
            sigma.deriv[y.dim] <- -dsigmapq %*% sigmaqq.inv %*% mus[x.dim]
            for (i in 1:n.x){
                dsigmapq <- numeric(n.x)
                dsigmapq[i] <- get.corpar(rhos, y.dim, x.dim[i])*sigmas[y.dim]
                dsigmaqq <- matrix(0, nrow = n.x, ncol = n.x)
                dsigmaqq[i, i] <- 2*sigmas[x.dim[i]]
                for (j in (1:n.x)[-i]){
                    dsigmaqq[i, j] <- dsigmaqq[j, i] <- get.corpar(rhos, x.dim[i], x.dim[j])*sigmas[x.dim[j]]
                }
                dsigmaqqinv <- -sigmaqq.inv %*% dsigmaqq %*% sigmaqq.inv
                sigma.deriv[x.dim[i]] <- -dsigmapq %*% sigmaqq.inv %*% mus[x.dim] -
                    sub.sigma.mat.y %*% dsigmaqqinv %*% mus[x.dim]
            }
            ## Now the rhos.
            rho.deriv <- numeric(length(rhos))
            for (i in 1:length(rhos)){
                dims <- as.numeric(strsplit(substr(names(rhos)[i], 4, nchar(names(rhos))), ",")[[1]])
                ## One dimension is for y.dim, the other is for an x.dim.
                if (any(dims == y.dim) & any(x.dim %in% dims)){
                    xx.dim <- dims[dims != y.dim]
                    dsigmapq <- numeric(n.x)
                    dsigmapq[xx.dim == x.dim] <- sigmas[y.dim]*sigmas[xx.dim]
                    rho.deriv[i] <- -dsigmapq %*% sigmaqq.inv %*% mus[x.dim]
                } else if (all(dims %in% x.dim)){
                    ## If both dimensions are for x.dims.
                    which.x.dim <- which(x.dim %in% dims)
                    dsigmaqq <- matrix(0, nrow = n.x, ncol = n.x)
                    dsigmaqq[which.x.dim[1], which.x.dim[2]] <-
                        dsigmaqq[which.x.dim[2], which.x.dim[1]] <- sigmas[dims[1]]*sigmas[dims[2]]
                    dsigmaqqinv <- -sigmaqq.inv %*% dsigmaqq %*% sigmaqq.inv
                    rho.deriv[i] <- -sub.sigma.mat.y %*% dsigmaqqinv %*% mus[x.dim]
                }
                ## Note that if a dimension is not an x.dim or a y.dim, the
                ## derivative is zero.
            }
            ## Filling in the matrix.
            betas.gr[1, ] <- c(mu.deriv, sigma.deriv, rho.deriv, rep(0, length(psis)), rep(0, length(phis)))

            ## Next we need the partial derivatives for the
            ## coefficients for the predictor variables, starting with
            ## respect to mu, which are just 0.
            mu.deriv <- matrix(0, nrow = n.x, ncol = length(mus))
            ## Now the sigmas.
            sigma.deriv <- matrix(0, nrow = n.x, ncol = length(sigmas))
            dsigmapq <- matrix(0, nrow = 1, ncol = n.x)
            for (i in 1:n.x){
                dsigmapq[1, i] <- get.corpar(rhos, y.dim, x.dim[i])*sigmas[x.dim[i]]
            }
            sigma.deriv[, y.dim] <- dsigmapq %*% sigmaqq.inv
            for (i in 1:n.x){
                dsigmapq <- numeric(n.x)
                dsigmapq[i] <- get.corpar(rhos, y.dim, x.dim[i])*sigmas[y.dim]
                dsigmaqq <- matrix(0, nrow = n.x, ncol = n.x)
                dsigmaqq[i, i] <- 2*sigmas[x.dim[i]]
                for (j in (1:n.x)[-i]){
                    dsigmaqq[i, j] <- dsigmaqq[j, i] <- get.corpar(rhos, x.dim[i], x.dim[j])*sigmas[x.dim[j]]
                }
                dsigmaqqinv <- -sigmaqq.inv %*% dsigmaqq %*% sigmaqq.inv
                sigma.deriv[, x.dim[i]] <- dsigmapq %*% sigmaqq.inv +
                    sub.sigma.mat.y %*% dsigmaqqinv
            }
            ## Now with respect to rho.
            rho.deriv <- matrix(0, nrow = n.x, ncol = length(rhos))
            for (i in 1:length(rhos)){
                dims <- as.numeric(strsplit(substr(names(rhos)[i], 4, nchar(names(rhos))), ",")[[1]])
                ## One dimension is for y.dim, the other is for an x.dim.
                if (any(dims == y.dim) & any(x.dim %in% dims)){
                    xx.dim <- dims[dims != y.dim]
                    dsigmapq <- numeric(n.x)
                    dsigmapq[xx.dim == x.dim] <- sigmas[y.dim]*sigmas[xx.dim]
                    rho.deriv[, i] <- dsigmapq %*% sigmaqq.inv
                } else if (all(dims %in% x.dim)){
                    ## If both dimensions are for x.dims.
                    which.x.dim <- which(x.dim %in% dims)
                    dsigmaqq <- matrix(0, nrow = n.x, ncol = n.x)
                    dsigmaqq[which.x.dim[1], which.x.dim[2]] <-
                        dsigmaqq[which.x.dim[2], which.x.dim[1]] <- sigmas[dims[1]]*sigmas[dims[2]]
                    dsigmaqqinv <- -sigmaqq.inv %*% dsigmaqq %*% sigmaqq.inv
                    rho.deriv[, i] <- sub.sigma.mat.y %*% dsigmaqqinv 
                }
            }
            ## Putting all the partial derivatives together.
            psiphi.deriv <- matrix(0, nrow = n.x, ncol = length(psis) + length(phis))
            betas.gr[-1, ] <- cbind(mu.deriv, sigma.deriv, rho.deriv, psiphi.deriv)
        } else {
            betas.gr <- matrix(0, nrow = 2, ncol = length(par.vec))
            ## Partial derivatives for the intercept.
            betas.gr[1, names(est) == paste0("mu", y.dim)] <- 1
            betas.gr[1, names(est) == paste0("mu", x.dim)] <- -sigmas[y.dim]/sigmas[x.dim]
            betas.gr[1, names(est) == paste0("sigma", y.dim)] <- -mus[x.dim]/sigmas[x.dim]
            betas.gr[1, names(est) == paste0("sigma", x.dim)] <- sigmas[y.dim]*mus[x.dim]/(sigmas[x.dim]^2)
            ## Partial derivatives for the slope.
            betas.gr[2, names(est) == paste0("sigma", y.dim)] <- 1/sigmas[x.dim]
            betas.gr[2, names(est) == paste0("sigma", x.dim)] <- -sigmas[y.dim]/(sigmas[x.dim]^2)
        }
        pars.varcov <- fit$vcov$varcov
        betas.varcov <- betas.gr %*% pars.varcov %*% t(betas.gr)
        betas.se <- sqrt(diag(betas.varcov))
        if (vcov){
            out <- list(est = betas.est, varcov = betas.varcov)
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

#' Predict Measurements from True Values
#'
#' @description
#' Makes predictions for one dimension based on true measurements of other dimensions.
#' Can handle multiple predictions simultaneously.
#'
#' @param object A fitted model object from fit.morph()
#' @param y.dim Integer specifying which dimension to predict
#' @param newdata A data frame of other dimensions to use for prediction. Column names
#'        must be of the form "dimX" where X is the dimension number. Can contain
#'        multiple rows for multiple predictions.
#' @param type Either "lm" or "pca" for prediction type
#' @param ... Additional arguments passed to methods
#'
#' @return A matrix with Estimate and Std.Error columns for each prediction
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

  # Check for NA values
  if (any(sapply(newdata, function(x) any(is.na(x))))) {
    stop("newdata cannot contain missing values (NA)")
  }

  # Get predictor dims
  x.dim <- as.numeric(sapply(strsplit(colnames(newdata), "dim"), function(x) x[2]))

  # Get coeficients and vcov matrix
  betas <- calc.betas(fit = object, y.dim = y.dim, x.dim = x.dim, type = type, vcov = FALSE)
  vcov.obj <- calc.betas(fit = object, y.dim = y.dim, x.dim = x.dim, type = type, vcov = TRUE)
  vcov.mat <- vcov.obj$varcov

  # Create design matrix for all preds
  n_preds <- nrow(newdata)
  X <- cbind(1, as.matrix(newdata))

  # Calculate preds for all rows
  pred.est <- X %*% betas[,"Estimate"]

  # Calculate se for all predictions
  pred.se <- sqrt(diag(X %*% vcov.mat %*% t(X)))

  # Combine
  result <- cbind(pred.est, pred.se)
  colnames(result) <- c("Estimate", "Std. Error")

  return(result)
}

# -------------------------------------------------------------------------------------------------------

#' Predict from Observed Measurements
#'
#' @description
#' Makes predictions by combining prior information with observed measurements that contain error.
#'
#' @param object A fitted model object from fit.morph()
#' @param ... Additional arguments:
#'        \itemize{
#'          \item true: A vector of true dimension measurements. Set elements to NA for dimensions to predict
#'          \item obs: A matrix of observed measurements with error, where each row represents measurements
#'                from one photograph. NA values allowed for missing measurements.
#'        }
#'
#' @return A numeric vector of predictions for all dimensions
#' @export
predict.obs <- function(object, ...) {
  # Extract arguments from ...
  args <- list(...)

  # If first argument isn't named, assume it's 'true'
  if (length(args) > 0 && is.null(names(args)[1])) {
    true <- args[[1]]
  } else {
    true <- args$true
  }

  # If second argument isn't named, assume it's 'obs'
  if (length(args) > 1 && is.null(names(args)[2])) {
    obs <- args[[2]]
  } else {
    obs <- args$obs
  }

  # Input validation
  if (!inherits(object, "lme.morph")) {
    stop("Invalid model object. Input must be a fitted model of class 'lme.morph'")
  }

  if (is.null(true)) {
    stop("'true' argument must be provided")
  }

  # Rest of function remains the same...
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
        out <- out + dmvnorm(obs[j, obs.dims],
                             mean = t[obs.dims],
                             sigma = pars$xi[obs.dims, obs.dims, drop = FALSE],
                                      log = TRUE)
      }
    }

    # Add prior
    out <- out + dmvnorm(t, mean = pars$mus, sigma = pars$sigma, log = TRUE)
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

## A function to grab a correlation parameter from a named vector.
get.corpar <- function(cors, i, j){
    inds <- substr(names(cors), 4, nchar(names(cors)))
    str <- paste(sort(c(i, j)), collapse = ",")
    cors[inds == str]
}
