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
#' @importFrom nlme getData
#' @importFrom lmeInfo extract_varcomp Fisher_info
#' @importFrom msm deltamethod
#' @importFrom Matrix bdiag
#'
#' @keywords internal
vcov.lme.morph <- function(object) {
  # Extract data from fitted model
  data <- getData(object)

  # Num of dimes in the morphometric data
  m <- length(unique(data$dim))

  # Get parameter estimates + names
  beta <- object$coefficients$fixed
  names(beta) <- paste0("beta", seq_along(beta))
  vc <- extract_varcomp(object, vector = TRUE)

  # Set up parameter vectors
  orig.est <- c(beta, vc)
  param_names <- names(orig.est)

  # Create dimension pairs for correlation parameters (only for dimensions we have)
  pairs <- combn(m, 2)
  dim.pairs <- apply(pairs, 2, function(x) paste(x, collapse = ","))

  # Initialise est vector and g list with correct dimensions
  est <- numeric(m + m + ncol(pairs) + m + ncol(pairs))
  names(est) <- c(paste0("mu", 1:m),
                  paste0("sigma", 1:m),
                  paste0("rho", dim.pairs),
                  paste0("psi", 1:m),
                  paste0("phi", dim.pairs))
  g <- vector("list", length(est))
  names(g) <- names(est)

  # Transform parameters
  k <- 1

  # 1. Mean parameters (mu)
  for(i in 1:m) {
    est[k] <- beta[i]
    g[[k]] <- reformulate(paste0("x", i))
    k <- k + 1
  }

  # 2. Variance parameters (sigma from Tau)
  for(i in 1:m) {
    tau_name <- paste0("Tau.animal.id.animal.id.var(dim", i, ")")
    if(tau_name %in% names(orig.est)) {
      tau_idx <- which(param_names == tau_name)
      est[k] <- sqrt(orig.est[tau_idx])
      g[[k]] <- reformulate(paste0("sqrt(x", tau_idx, ")"))
    } else {
      est[k] <- NA
      g[[k]] <- reformulate("0")
    }
    k <- k + 1
  }

  # 3. Corrrelation parameters (rho)
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      var_i <- paste0("Tau.animal.id.animal.id.var(dim", i, ")")
      var_j <- paste0("Tau.animal.id.animal.id.var(dim", j, ")")
      cov_ij <- paste0("Tau.animal.id.animal.id.cov(dim", j, ",dim", i, ")")

      if(all(c(var_i, var_j, cov_ij) %in% names(orig.est))) {
        var_i_idx <- which(param_names == var_i)
        var_j_idx <- which(param_names == var_j)
        cov_ij_idx <- which(param_names == cov_ij)

        est[k] <- orig.est[cov_ij_idx] /
          (sqrt(orig.est[var_i_idx]) * sqrt(orig.est[var_j_idx]))

        g[[k]] <- reformulate(paste0("x", cov_ij_idx, "/(sqrt(x", var_i_idx,
                                     ")*sqrt(x", var_j_idx, "))"))
      } else {
        est[k] <- NA
        g[[k]] <- reformulate("0")
      }
      k <- k + 1
    }
  }

  # 4. Standard deviation parameters (psi)
  sigma_sq_name <- "sigma_sq"
  if(sigma_sq_name %in% names(orig.est)) {
    sigma_sq_idx <- which(param_names == sigma_sq_name)
    for(i in 1:m) {
      if(i == 1) {
        est[k] <- sqrt(orig.est[sigma_sq_idx])
        g[[k]] <- reformulate(paste0("sqrt(x", sigma_sq_idx, ")"))
      } else {
        var_param <- paste0("var_params", i-1)
        if(var_param %in% names(orig.est)) {
          var_idx <- which(param_names == var_param)
          est[k] <- sqrt(orig.est[sigma_sq_idx]) * orig.est[var_idx]
          g[[k]] <- reformulate(paste0("sqrt(x", sigma_sq_idx, ")*x", var_idx))
        } else {
          est[k] <- NA
          g[[k]] <- reformulate("0")
        }
      }
      k <- k + 1
    }
  } else {
    for(i in 1:m) {
      est[k] <- NA
      g[[k]] <- reformulate("0")
      k <- k + 1
    }
  }

  # 5. Correlation parameters (phi)
  l <- 1
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      cor_param <- paste0("cor_params", l)
      if(cor_param %in% names(orig.est)) {
        cor_idx <- which(param_names == cor_param)
        est[k] <- orig.est[cor_idx]
        g[[k]] <- reformulate(paste0("x", cor_idx))
      } else {
        est[k] <- NA
        g[[k]] <- reformulate("0")
      }
      l <- l + 1
      k <- k + 1
    }
  }

  # Convert to standard lme class for vcov extraction
  object.lme <- object
  class(object.lme) <- "lme"

  # Get original ariance-covariance matrix using Fisher info
  orig.varcov <- as.matrix(bdiag(vcov(object.lme),
                                 solve(Fisher_info(object))))

  # Transform variance-covariance matrix using delta metod
  varcov <- deltamethod(g, orig.est, orig.varcov, ses = FALSE)
  rownames(varcov) <- colnames(varcov) <- names(est)

  # Return estimates and variance-covariance matrix
  list(est = est, varcov = varcov)
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










