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
#' @param x A data frame containing:
#'   \itemize{
#'     \item animal.id: Factor indicating the animal
#'     \item photo.id: Factor indicating the photo
#'     \item dim: Factor indicating the dimension
#'     \item measurement: The observed measurement
#'   }
#' @param ... Additional arguments:
#'   \itemize{
#'     \item dims: Integer vector of length 2, selecting which dimensions to plot
#'     \item plot.data: Logical; if FALSE, the plotting area is set up but points aren't plotted
#'     \item xlim,ylim: Numeric vectors of length 2 giving plot limits
#'     \item ratios: Logical; if TRUE, y-axis represents the ratio between dimensions
#'     \item xlab,ylab: Character strings giving axis labels
#'   }
#'
#' @return NULL (invisibly). Creates a plot as a side effect.
#'
#' @export
plot.morph <- function(x, ...) {
  # Extract arguments from ...
  args <- list(...)
  dims <- if (!is.null(args$dims)) args$dims else c(1, 2)
  plot.data <- if (!is.null(args$plot.data)) args$plot.data else TRUE
  xlim <- args$xlim
  ylim <- args$ylim
  ratios <- if (!is.null(args$ratios)) args$ratios else FALSE
  xlab <- args$xlab
  ylab <- args$ylab

  # Input validation
  if (!is.data.frame(x)) {
    stop("Columns missing from data frame. See ?plot.morph for column requirements.")
  }

  # Check required columns exist
  required_cols <- c("animal.id", "photo.id", "dim", "measurement")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in data frame:\n",
      paste0("  - '", missing_cols, "'", collapse = "\n")
    )
  }

  # Check dims
  if (!is.numeric(dims) || length(dims) != 2) {
    stop("'dims' must be a numeric vector of length 2")
  }

  # Check if all requested dims exist in data
  if (!all(dims %in% as.numeric(levels(x$dim)))) {
    stop("Not all requested dimensions are present in the data")
  }

  # Keep only the dims to be plotted
  x <- x[x$dim %in% dims, ]

  # Check if we have both dimensions
  if (!all(dims %in% unique(x$dim))) {
    stop("Not all requested dimensions are present in the data")
  }

  # Extract columns
  measurement <- x$measurement
  dim <- x$dim
  animal.id <- x$animal.id
  photo.id <- x$photo.id

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
#' @export
plot.lme.morph <- function(x, dims = c(1, 2), type = "data",
                           line.type = "none", confints = !add,
                           add = FALSE, reverse.axes = FALSE,
                           plot.data = TRUE, xlim = NULL, ylim = NULL,
                           xlab = NULL, ylab = NULL, ...) {
  # Valid types
  valid_types <- c("data", "ratio")
  if (!type %in% valid_types) {
    stop(
      "Invalid type argument. See ?plot.lme.morph for possible selections."
    )
  }

  # Valid line types
  valid_line_types <- c("none", "lm", "pca")
  if (!line.type %in% valid_line_types) {
    stop(
      "Invalid line.type argument. See ?plot.lme.morph for possible selections."
    )
  }

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
            predict(x, y.dim = dims[2], true_measurements = newdata)
          } else {
            predict(x, y.dim = dims[2], true_measurements = newdata, type = "pca")
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
#'   pars = c(315, 150, 100,    # means
#'            25, 15, 10,       # SDs
#'            0.85, 0.80, 0.75,  # correlations
#'            10, 6, 4,    # measurement SDs
#'            0.5, 0.4, 0.3)  # measurement correlations
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
    stop(
      "Invalid model object. Input must be a fitted model of class 'lme.morph'.\n",
      "Use fit.morph() to create a valid model object."
    )
  }

  # Extract fixed effects for means
  mus <- nlme::fixef(x)

  # Extract variance components from random affects
  rand_effects <- x$modelStruct$reStruct[[1]]
  sigmas <- sqrt(diag(as.matrix(rand_effects)))
  # Get corelation from first off-diagonal element
  rhos <- as.matrix(rand_effects)[1,2] / (sigmas[1] * sigmas[2])

  # Check dims exist and provide informative error
  m <- length(mus)
  available_dims <- 1:m
  invalid_dims <- c()

  if (dim1 > m || dim1 < 1) invalid_dims <- c(invalid_dims, dim1)
  if (dim2 > m || dim2 < 1) invalid_dims <- c(invalid_dims, dim2)

  if (length(invalid_dims) > 0) {
    stop(
      "Invalid dimension(s) requested. Available dimensions are:\n",
      paste0("  - ", available_dims, collapse = "\n"),
      "\nYou requested: ", paste(invalid_dims, collapse = ", ")
    )
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
    stop(
      "Error computing ratio distribution.\n",
      "This might occur if:\n",
      "  - The dimensions have very different scales\n",
      "  - There is extremely high correlation\n",
      "  - The denominator dimension has values near zero\n"
    )
  })

  invisible(NULL)
}
