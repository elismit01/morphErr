#' Visualisation Functions for morphErr
#'
#' This file contains functions for visualising morphometric data and model results:
#' - plot.morph(): Plots raw morphometric data
#' - plot.lme.morph(): Creates model-based plots
#' - plot.ratio.pdf(): Plots probability density functions for ratios
#'
#' @name plot
#' @keywords internal
#' @importFrom graphics plot.new plot.window box axis title points lines par
#' @importFrom grDevices hcl.colors
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
#'   pars = c(100, 50, 25,    # means
#'            10, 5, 2,       # SDs
#'            0.7, 0.6, 0.5,  # correlations
#'            1, 0.5, 0.2,    # measurement SDs
#'            0.3, 0.2, 0.1)  # measurement correlations
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
#'
#' @examples
#' \dontrun{
#' # Generate some data
#' data <- sim.measurements(
#'   n.animals = 10,
#'   n.photos = rep(3, 10),
#'   m = 2,
#'   pars = c(100, 50, # means
#'            10, 5, # Sds
#'            0.7,# correlation
#'            1, 0.5,  # measurement SDs
#'            0.3) # measurement correlation
#' )
#'
#' # Fit model and create plots
#' fit <- fit.morph(data)
#' plot(fit, type = "data", line.type = "lm")
#' plot(fit, type = "ratio")
#' }
#'
#' @export
plot.lme.morph <- function(x, dims = c(1, 2), type = "data",
                           line.type = "none", confints = !add,
                           add = FALSE, reverse.axes = FALSE,
                           plot.data = TRUE, xlim = NULL, ylim = NULL,
                           xlab = NULL, ylab = NULL, ...){
  # Get data
  data <- nlme::getData(x)

  # Basic data plot if requested
  if (type == "data"){
    if (!add){
      plot.morph(data, dims, plot.data = plot.data,
                 xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab)
    }

    # Get plot limits if not provided
    xlim <- par("usr")[c(1, 2)]

    # Add lines if requested
    if (line.type != "none"){
      # Get coefficients + check they exist
      tryCatch({
        if (line.type == "lm"){
          betas <- summary(x, type = "betas",
                           y.dim = dims[2], x.dim = dims[1])
        } else if (line.type == "pca"){
          betas <- summary(x, type = "betas-pca",
                           y.dim = dims[2], x.dim = dims[1])
        }

        if (!is.matrix(betas) || nrow(betas) < 2) {
          stop("Could not compute valid coefficients")
        }

        # Draw line
        if (reverse.axes){
          abline(-betas[1, 1]/betas[2, 1], 1/betas[2, 1], ...)
        } else {
          abline(betas[1:2, 1], ...)
        }

        # Add confidence intervals if request
        if (confints){
          dim1.xx <- seq(xlim[1], xlim[2], length.out = 1000)
          newdata <- data.frame(dim1.xx)
          colnames(newdata) <- paste0("dim", dims[1])

          if (line.type == "lm"){
            preds <- predict(x, y.dim = dims[2], newdata = newdata)
          } else {
            preds <- predict(x, y.dim = dims[2],
                             newdata = newdata, type = "pca")
          }

          preds.upper <- preds[, 1] + qnorm(0.975)*preds[, 2]
          preds.lower <- preds[, 1] - qnorm(0.975)*preds[, 2]
          lines(dim1.xx, preds.upper, lty = "dotted")
          lines(dim1.xx, preds.lower, lty = "dotted")
        }
      }, error = function(e) {
        warning("Could not add line to plot: ", e$message)
      })
    }
  } else if (type == "ratio"){
    if (!add){
      plot.morph(data, dims, xlim = xlim, ylim = ylim,
                 ratios = TRUE, plot.data = plot.data,
                 xlab = xlab, ylab = ylab)
    }

    # plot limits if not specified
    xlim <- par("usr")[c(1, 2)]

    if (line.type != "none"){
      xx <- seq(xlim[1], xlim[2], length.out = 1000)
      preds <- calc.conditional.ratio(x, y.dim = dims[2],
                                      x.dim = dims[1],
                                      newdata.x.dim = xx,
                                      type = line.type)
      lines(xx, preds[, 1])

      if (confints){
        lines(xx, preds[, 1] + qnorm(0.975)*preds[, 2],
              lty = "dotted")
        lines(xx, preds[, 1] - qnorm(0.975)*preds[, 2],
              lty = "dotted")
      }
    }
  }

  invisible(NULL)
}

