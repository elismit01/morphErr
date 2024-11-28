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

