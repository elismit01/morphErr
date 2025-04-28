# Visualisation Functions for morphErr
#
# This file contains functions for visualising morphometric data and model results:
# - plotmorph(): Plots raw morphometric data
# - plot.lme.morph(): Creates model-based plots
# - plot.ratio.pdf(): Plots probability density functions for ratios

#' Plot Morphometric Data
#'
#' Creates scatter plots of morphometric measurements, with options for
#' different dimension combinations and plotting of ratios.
#'
#' @section The `data` argument:
#' 
#'  The [`manta`] object is an example of a correctly formatted
#' \code{data} argument. It must be a data frame with the following
#' columns:
#'
#' \describe{
#' 
#'  \item{\code{animal.id}}{An individual identification number. Rows
#'                    with the same \code{animal.id} correspond to
#'                    measurements of the same individual manta ray.}
#'
#'  \item{\code{photo.id}}{A photo identification number. Rows with
#'                   the same \code{photo.id} correspond to
#'                   measurements taken from the same image.}
#'
#'  \item{\code{photo.id}}{An integer indicating the dimension the
#'              measurement is for.}
#'
#'  \item{\code{measurement}}{The corresponding measurement.}
#' 
#' }
#'
#' @param data A data frame containing the morphometric data. See the
#'     section below on the correct formatting of this argument.
#' @param dims Integer vector of length 2, selecting which dimensions
#'     to plot.
#' @param plot.data Logical. If `FALSE`, the plotting area is set up
#'     but points aren't plotted.
#' @param ratios Logical. If `TRUE`, the y-axis represents the ratio
#'     between dimensions.
#' @param xlim,ylim Limits for the axes.
#' @param xlab,ylab Titles for the axes.
#'
#' @examples
#' ## Plotting dimensions 1 and 2.
#' plotmorph(manta)
#' ## Plotting dimensions 1 and 3.
#' plotmorph(manta, dims = c(1, 3))
#' ## Plotting the ratio of dimension 2 divided by dimension 1 on the
#' ## y-axis.
#' plotmorph(manta, dims = c(1, 2), ratios = TRUE)
#' 
#' @export
plotmorph <- function(data, dims = c(1, 2), plot.data = TRUE, ratios = FALSE,
                      xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL){

  # Input validation
  if (!is.data.frame(data)) {
    stop("The data argument must be a data frame. See ?plotmorph for column requirements.")
  }

  # Turning dim into a factor if it is isn't already.
  if (!is.factor(data$dim)){
    data$dim <- factor(data$dim)
  }

  # Check required columns exist
  required_cols <- c("animal.id", "photo.id", "dim", "measurement")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in data:\n",
      paste0("  - '", missing_cols, "'", collapse = "\n")
    )
  }

  # Check dims
  if (!is.numeric(dims) || length(dims) != 2) {
    stop("'dims' must be a numeric vector of length 2")
  }

  # Check if all requested dims exist in data
  if (!all(dims %in% as.numeric(levels(data$dim)))) {
    stop("Not all requested dimensions are present in the data")
  }

  # Keep only the dims to be plotted
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
  cols <- hcl.colors(n.animals, palette = "Dark 2")

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
  data <- getData(x)

  # Basic data plot if requested
  if (type == "data") {
    if (!add) {
      plotmorph(data, plot.data = plot.data,
                 xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab,
                 dims = dims)
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
            predict(x, y.dim = dims[2], newdata = newdata)
          } else {
            predict(x, y.dim = dims[2], newdata = newdata, type = "pca")
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
      plotmorph(data, xlim = xlim, ylim = ylim,
                 ratios = TRUE, plot.data = plot.data,
                 xlab = xlab, ylab = ylab,
                 dims = dims)
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
