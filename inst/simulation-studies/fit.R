# Load morphErr
library(morphErr)

#' Helper Functions for Model Fitting
#'
#' This script contains functions for fitting models to both
#' simulated data and (eventually) real manta ray data.

#' Create example fit object with simulated data
create_example_fit <- function() {
  # Simulate data using params from paper
  n.animals <- 50
  n.photos <- 5

  cat("Creating data with:", n.animals, "animals and", n.photos, "photos each\n")

  # Create param vector
  pars <- c(290, 125, 75,    # means
            45, 25, 15,       # SDs
            0.75, 0.80, 0.85, # correlations
            2.0, 1.5, 1.0,    # measurement SDs
            0.4, 0.5, 0.6)    # measurement correlations

  # Print true parames
  cat("\nTrue Parameters:\n")
  cat("Means:", pars[1:3], "\n")
  cat("SDs:", pars[4:6], "\n")
  cat("Correlations:", pars[7:9], "\n")
  cat("Scale params:", pars[10:12], "\n")
  cat("Error params:", pars[13:15], "\n\n")

  test_data <- sim.measurements(
    n.animals = n.animals,
    n.photos = rep(n.photos, n.animals),
    m = 3,
    pars = pars
  )

  # Check data structure
  cat("Data structure:\n")
  print(str(test_data))

  # Fit model
  cat("\nFitting model...\n")
  fit <- suppressWarnings(fit.morph(test_data))

  return(fit)
}

#' Fit modle to manta ray data
fit_manta_model <- function(data_path = NULL) {
  # If no data provided, use example data
  if (is.null(data_path)) {
    fit <- create_example_fit()
  } else {
    # TODO: Add code to load and process real manta ray data
    stop("Real data loading not yet implemented")
  }

  # Save
  save(fit, file = "manta-fit.RData")
  return(fit)
}

# Run if script is executed directly
if (sys.nframe() == 0) {
  # Create example fit object
  cat("Creating new fit object...\n")
  fit <- create_example_fit()

  # Print summaries
  cat("\nModel Summary:\n")
  print(summary(fit))

  cat("\nIsometry Tests:\n")
  print(summary(fit, type = "isometric-pca"))

  cat("\nExample fit object 'fit' created in environment\n")
}
