#' Visualisation Functions for morphErr
#'
#' This file contains functions for visualising morphometric data and model results:
#' - plot.morph(): Plots raw morphometric data
#' - plot.lme.morph(): Creates model-based plots
#' - plot.ratio.pdf(): Plots probability density functions for ratios
#'
#' @name plot
#' @keywords internal
NULL

#' Plot Morphometric Data
#'
#' Creates scatter plots of morphometric measurements, with options for
#' different dimension combinations and ratio calculations.
#'
#' @param data A data frame containing:
#'   \itemize{
#'     \item animal.id: Factor indicating the animal
#'     \item photo.id: Factor indicating the photo
#'     \item dim: Factor indicating the dimension
#'     \item measurement: The observed measurement
#'   }
#' @param dims Integer vector of length 2, selecting which dimensions to plot
#' @param plot.data Logical; if FALSE, the plotting area is set up but points aren't plotted
#' @param xlim,ylim Numeric vectors of length 2 giving plot limits
#' @param ratios Logical; if TRUE, y-axis represents the ratio between dimensions
#' @param xlab,ylab Character strings giving axis labels
#'
#' @return NULL (invisibly). Creates a plot as a side effect.
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' data <- sim.measurements(
#'   n.animals = 5,
#'   n.photos = rep(3, 5),
#'   m = 3,
#'   pars = c(290, 125, 75,    # means
#'            45, 25, 15,       # SDs
#'            0.75, 0.80, 0.85,  # correlations
#'            2.0, 1.5, 1.0,    # measurement SDs
#'            0.4, 0.5, 0.6)  # measurement correlations
#' )
#'
#' # Basic scatter plot
#' plot.morph(data, dims = c(1, 2))
#'
#' # Plot with ratios
#' plot.morph(data, dims = c(1, 2), ratios = TRUE)
#' }
#'
#' @export
plot.morph <- function(data, dims = c(1, 2), plot.data = TRUE,
                       xlim = NULL, ylim = NULL,
                       ratios = FALSE, xlab = NULL, ylab = NULL){
  # Keep only the dimensions to be plotted
  data <- data[data$dim %in% dims, ]

  # Check if we have both dimensions
  if (!all(dims %in% unique(data$dim))) {
    stop("Not all requested dimensions are present in the data")
  }

  # Extract columns
  measurement <- data$measurement
  dim <- data$dim
  animal.id <- data$animal.id
  photo.id <- data$photo.id

  # Total number of animals
  n.animals <- length(unique(animal.id))

  # Choose colour palette
  cols <- hcl.colors(n.animals, palette = "viridis")

  # Map colours to ids
  col.id.map <- data.frame(cols, unique(animal.id))

  # Create plotting area
  plot.new()

  # Calculate plot limits if not provided
  if (is.null(ylim)){
    if (ratios){
      ylim <- range(measurement[dim == dims[2]]/measurement[dim == dims[1]])
    } else {
      ylim <- range(measurement[dim == dims[2]])
    }
  }
  if (is.null(xlim)){
    xlim <- range(measurement[dim == dims[1]])
  }

  # Set up plot window
  plot.window(xlim = xlim, ylim = ylim)

  # Add box and axes
  box()
  axis(1)
  axis(2)

  # Set up labells
  if (ratios){
    if (is.null(xlab)) xlab <- paste0("dim", dims[1])
    if (is.null(ylab)) ylab <- paste0("dim", dims[2], "/dim", dims[1])
  } else {
    if (is.null(xlab)) xlab <- paste0("dim", dims[1])
    if (is.null(ylab)) ylab <- paste0("dim", dims[2])
  }
  title(xlab = xlab, ylab = ylab)

  # Plot points if requested
  if (plot.data){
    for (i in unique(animal.id)){
      photo.ids <- unique(photo.id[animal.id == i])
      col <- col.id.map[col.id.map[, 2] == i, 1]
      for (j in photo.ids){
        if (ratios){
          points(measurement[dim == dims[1] & animal.id == i &
                               photo.id == j],
                 measurement[dim == dims[2] & animal.id == i &
                               photo.id == j]/
                   measurement[dim == dims[1] & animal.id == i &
                                 photo.id == j],
                 col = col, pch = 16)
        } else {
          points(measurement[dim == dims[1] & animal.id == i &
                               photo.id == j],
                 measurement[dim == dims[2] & animal.id == i &
                               photo.id == j],
                 col = col, pch = 16)
        }
      }
    }
  }

  invisible(NULL)
}

# -------------------------------------------------------------------------------------------------------

#' Plot Method for Morphometric Model Fits
#'
#' S3 method for plotting fitted morphometric models. Supports various plot types
#' and can overlay fitted lines on data.
#'
#' @param x An object of class "lme.morph" returned by fit.morph()
#' @param dims Integer vector of length 2, selecting which dimensions to plot
#' @param type Character string specifying plot type:
#'   \itemize{
#'     \item "data": plots the raw data
#'     \item "ratio": plots ratios between dimensions
#'   }
#' @param line.type Character string specifying type of line to overlay:
#'   \itemize{
#'     \item "none": no line
#'     \item "lm": linear regression line
#'     \item "pca": principal components line
#'   }
#' @param confints Logical; if TRUE, includes confidence intervals
#' @param add Logical; if TRUE, adds to existing plot
#' @param reverse.axes Logical; if TRUE, swaps x and y interpretations
#' @param plot.data Logical; if FALSE, only plots fitted lines
#' @param xlim,ylim Numeric vectors of length 2 giving plot limits
#' @param xlab,ylab Character strings giving axis labels
#' @param ... Additional arguments passed to plotting functions
#'
#' @return NULL (invisibly). Creates a plot as a side effect.
plot.lme.morph <- function(x, dims = c(1, 2), type = "data",
                           line.type = "none", confints = !add,
                           add = FALSE, reverse.axes = FALSE,
                           plot.data = TRUE, xlim = NULL, ylim = NULL,
                           xlab = NULL, ylab = NULL, ...) {
  # Get data
  data <- nlme::getData(x)

  # Basic data plot if requested
  if (type == "data") {
    if (!add) {
      plot.morph(data, dims, plot.data = plot.data,
                 xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab)
    }

    # Get plot limits if not provided
    xlim <- par("usr")[c(1, 2)]

    # Add lines if requested
    if (line.type != "none") {
      # Get coefficients + check they exist
      suppressWarnings({
        if (line.type == "lm") {
          betas <- summary(x, type = "betas",
                           y.dim = dims[2], x.dim = dims[1])
        } else if (line.type == "pca") {
          betas <- summary(x, type = "betas-pca",
                           y.dim = dims[2], x.dim = dims[1])
        }

        if (!is.matrix(betas) || nrow(betas) < 2) {
          return(invisible())
        }

        # Draw line (if coefs valid)
        if (!is.na(betas[1,1]) && !is.na(betas[2,1]) && betas[2,1] != 0) {
          # First main line
          if (reverse.axes) {
            abline(-betas[1, 1]/betas[2, 1], 1/betas[2, 1], ...)
          } else {
            abline(betas[1:2, 1], ...)
          }

          # Add confidence intervals if request
          if (confints) {
            # Sequence of x values
            dim1.xx <- seq(xlim[1], xlim[2], length.out = 1000)

            # Prediction data frame
            newdata <- data.frame(x = dim1.xx)
            names(newdata) <- paste0("dim", dims[1])

            # Get preds
            preds <- if (line.type == "lm") {
              suppressWarnings(predict(x, y.dim = dims[2], true_measurements = newdata))
            } else {
              suppressWarnings(predict(x, y.dim = dims[2], true_measurements = newdata, type = "pca"))
            }

            # Only add lines if gt valid preds
            if (is.matrix(preds) && nrow(preds) == length(dim1.xx)) {
              preds.upper <- preds[, 1] + qnorm(0.975)*preds[, 2]
              preds.lower <- preds[, 1] - qnorm(0.975)*preds[, 2]
              lines(dim1.xx, preds.upper, lty = "dotted")
              lines(dim1.xx, preds.lower, lty = "dotted")
            }
          }
        }
      })
    }
  } else if (type == "ratio") {
    if (!add) {
      plot.morph(data, dims, xlim = xlim, ylim = ylim,
                 ratios = TRUE, plot.data = plot.data,
                 xlab = xlab, ylab = ylab)
    }

    # plot limits if not specified
    xlim <- par("usr")[c(1, 2)]

    if (line.type != "none") {
      xx <- seq(xlim[1], xlim[2], length.out = 1000)
      suppressWarnings({
        preds <- calc.conditional.ratio(x, y.dim = dims[2],
                                        x.dim = dims[1],
                                        newdata.x.dim = xx,
                                        type = line.type)
        if (is.matrix(preds)) {
          lines(xx, preds[, 1])

          if (confints) {
            lines(xx, preds[, 1] + qnorm(0.975)*preds[, 2],
                  lty = "dotted")
            lines(xx, preds[, 1] - qnorm(0.975)*preds[, 2],
                  lty = "dotted")
          }
        }
      })
    }
  }

  invisible(NULL)
}

# -------------------------------------------------------------------------------------------------------

#' Plot Probability Density Function for Morphometric Ratios
#'
#' Creates a plot of the probability density function for the ratio
#' between two dimensions.
#'
#' @param x An object of class "lme.morph" returned by fit.morph()
#' @param dim1,dim2 Integers specifying which dimensions to use
#'
#' @return NULL (invisibly). Creates a plot as a side effect.
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' data <- sim.measurements(
#'   n.animals = 10,
#'   n.photos = rep(3, 10),
#'   m = 3,
#'   pars = c(290, 125, 75,    # means
#'            45, 25, 15,       # SDs
#'            0.75, 0.80, 0.85,  # correlations
#'            2.0, 1.5, 1.0,    # measurement SDs
#'            0.4, 0.5, 0.6)  # measurement correlations
#' )
#'
#' # Fit model and plot ratio PDF
#' fit <- fit.morph(data)
#' plot.ratio.pdf(fit, 1, 2)
#' }
#'
#' @keywords internal
plot.ratio.pdf <- function(x, dim1, dim2) {
  # Check if x is a valid model fit
  if (!inherits(x, "lme.morph")) {
    stop("'x' must be an object of class 'lme.morph'")
  }

  # Extract fixed effects for means
  mus <- nlme::fixef(x)

  # Extract variance components from random affects
  rand_effects <- x$modelStruct$reStruct[[1]]
  sigmas <- sqrt(diag(as.matrix(rand_effects)))
  # Get corelation from first off-diagonal element
  rhos <- as.matrix(rand_effects)[1,2] / (sigmas[1] * sigmas[2])

  # Check dimensions exist
  m <- length(mus)
  if (dim1 > m || dim2 > m) {
    stop("Requested dimensions exceed available dimensions")
  }

  # Construct correlation matrix
  sigma.mat <- diag(sigmas^2)
  sigma.mat[1, 2] <- sigma.mat[2, 1] <- rhos*sigmas[1]*sigmas[2]
  cor.mat <- cov2cor(sigma.mat)

  # calculating ratio distribution parameters
  tryCatch({
    mean.ratio <- mean.rcnorm.approx(mus[dim1], mus[dim2],
                                     sqrt(sigma.mat[dim1, dim1]),
                                     sqrt(sigma.mat[dim2, dim2]),
                                     cor.mat[dim1, dim2])
    sd.ratio <- sd.rcnorm.approx(mus[dim1], mus[dim2],
                                 sqrt(sigma.mat[dim1, dim1]),
                                 sqrt(sigma.mat[dim2, dim2]),
                                 cor.mat[dim1, dim2])

    # Create plot
    xx <- seq(mean.ratio - 5*sd.ratio, mean.ratio + 5*sd.ratio,
              length.out = 1000)
    yy <- drcnorm.approx(xx, mus[dim1], mus[dim2],
                         sqrt(sigma.mat[dim1, dim1]),
                         sqrt(sigma.mat[dim2, dim2]),
                         cor.mat[dim1, dim2])

    plot(xx, yy, type = "l",
         xlab = paste0("dim", dim1, "/dim", dim2),
         ylab = "Probability density")

  }, error = function(e) {
    stop("Could not compute ratio distribution: ", e$message)
  })

  invisible(NULL)
}

# -------------------------------------------------------------------------------------------------------

#' Plot Probability Density Function for Morphometric Ratios
#'
#' Creates a plot of the probability density function for the ratio
#' between two dimensions.
#'
#' @param x An object of class "lme.morph" returned by fit.morph()
#' @param dim1,dim2 Integers specifying which dimensions to use
#'
#' @return NULL (invisibly). Creates a plot as a side effect.
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' data <- sim.measurements(
#'   n.animals = 10,
#'   n.photos = rep(3, 10),
#'   m = 3,
#'   pars = c(290, 125, 75,    # means
#'            45, 25, 15,       # SDs
#'            0.75, 0.80, 0.85,  # correlations
#'            2.0, 1.5, 1.0,    # measurement SDs
#'            0.4, 0.5, 0.6)  # measurement correlations
#' )
#'
#' # Fit model and plot ratio PDF
#' fit <- fit.morph(data)
#' plot.ratio.pdf(fit, 1, 2)
#' }
#'
#' @keywords internal
plot.ratio.pdf <- function(x, dim1, dim2) {
  # Check if x is a valid model fit
  if (!inherits(x, "lme.morph")) {
    stop("'x' must be an object of class 'lme.morph'")
  }

  # Extract fixed effects for means
  mus <- nlme::fixef(x)

  # Extract variance components from random affects
  rand_effects <- x$modelStruct$reStruct[[1]]
  sigmas <- sqrt(diag(as.matrix(rand_effects)))
  # Get corelation from first off-diagonal element
  rhos <- as.matrix(rand_effects)[1,2] / (sigmas[1] * sigmas[2])

  # Check dimensions exist
  m <- length(mus)
  if (dim1 > m || dim2 > m) {
    stop("Requested dimensions exceed available dimensions")
  }

  # Construct correlation matrix
  sigma.mat <- diag(sigmas^2)
  sigma.mat[1, 2] <- sigma.mat[2, 1] <- rhos*sigmas[1]*sigmas[2]
  cor.mat <- cov2cor(sigma.mat)

  # calculating ratio distribution parameters
  tryCatch({
    mean.ratio <- mean.rcnorm.approx(mus[dim1], mus[dim2],
                                     sqrt(sigma.mat[dim1, dim1]),
                                     sqrt(sigma.mat[dim2, dim2]),
                                     cor.mat[dim1, dim2])
    sd.ratio <- sd.rcnorm.approx(mus[dim1], mus[dim2],
                                 sqrt(sigma.mat[dim1, dim1]),
                                 sqrt(sigma.mat[dim2, dim2]),
                                 cor.mat[dim1, dim2])

    # Create plot
    xx <- seq(mean.ratio - 5*sd.ratio, mean.ratio + 5*sd.ratio,
              length.out = 1000)
    yy <- drcnorm.approx(xx, mus[dim1], mus[dim2],
                         sqrt(sigma.mat[dim1, dim1]),
                         sqrt(sigma.mat[dim2, dim2]),
                         cor.mat[dim1, dim2])

    plot(xx, yy, type = "l",
         xlab = paste0("dim", dim1, "/dim", dim2),
         ylab = "Probability density")

  }, error = function(e) {
    stop("Could not compute ratio distribution: ", e$message)
  })

  invisible(NULL)
}


