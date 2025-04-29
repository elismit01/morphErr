# Simulation Functions for morphErr
#
# This file contains functions for simulating morphometric data:
# - sim.measurements(): Simulates measurement data from a survey
# - sim.morph(): Simulates multiple datasets and fits models
# - extract.sim.morph(): Extracts results from simulation studies


#' Simulate Morphometric Measurements
#'
#' Simulates morphometric data from a photogrammetry survey, with
#' observations subject to measurement error, under the model
#' described by \href{../doc/model.pdf}{Stevenson, Smit, and Setyawan
#' (in submission)}.
#'
#' @details For arguments `rhos` and `phis`, the elements must be
#'     ordered so that all \eqn{m - 1} correlations involving
#'     dimension 1 appear first in ascending numerical order, followed
#'     by all remaining \eqn{m - 2} correlations involving dimension
#'     2, and so on.
#'
#' For example, if `m = 4`, then the first three elements are the
#' correlations between dimension 1 and dimensions 2, 3, and 4,
#' respectively. The following two elements are correlations between
#' dimension 2 and dimensions 3 and 4, respectively. The final element
#' is the correlation between dimension 3 and 4.
#'
#' @param n.animals Integer. The number of animals in the sample.
#' @param n.photos Integer vector. If there are `n.animals` elements,
#'     each one specifies the number of photos for an individual. If
#'     there is one element, then that number of photos is used for
#'     all animals.
#' @param m Integer. The number of dimensions.
#' @param mus A vector with `m` elements providing the means of the
#'     true dimension sizes in the population.
#' @param sigmas A vector with `m` elements providing the standard
#'     deviations for true dimension sizes in the population.
#' @param rhos A vector, with one element for each pair of dimensions,
#'     providing the pairwise correlations between true dimension
#'     sizes in the population. See 'Details' for the correct order
#'     for the correlations.
#' @param psis A vector with `m` elements, providing the standard
#'     deviations of measurement errors for the dimensions.
#' @param phis A vector, with one element for each pair of dimensions,
#'     providing the pairwise correlations between measurement errors
#'     for the dimensions. See 'Details' for the correct order for the
#'     correlations.
#'
#' @seealso [`sim.morph()`] to conduct a simulation study by
#'     simulating multiple data sets and fitting a model to each one.
#'
#' @return
#' A data frame with four columns:
#' \describe{
#' 
#'   \item{\code{animal.id}}{An individual identification number. Rows
#'                    with the same \code{animal.id} correspond to
#'                    measurements of the same individual.}
#' 
#'   \item{\code{photo.id}}{A photo identification number. Rows with
#'                   the same \code{photo.id} correspond to
#'                   measurements taken from the same image.}
#' 
#'   \item{\code{dim}}{An integer indicating the dimension the
#'              measurement is for.}
#' 
#'   \item{\code{measurement}}{The observed measurement value.}
#' 
#' }
#'
#' @examples
#' ## Simulating data for ten animals, with two photos each, measuring
#' ## three dimensions.
#' data <- sim.measurements(n.animals = 10, n.photos = 2, m = 3,
#'                          mus = c(315, 150, 100),
#'                          sigmas = c(25, 15, 10),
#'                          rhos = c(0.85, 0.80, 0.75),
#'                          psis = c(10, 6, 4),
#'                          phis = c(0.5, 0.4, 0.3))
#' head(data)
#'
#' @export
sim.measurements <- function(n.animals, n.photos, m, mus, sigmas,
                             rhos, psis, phis){
  # If n.photos is scalar, repliicate for all animals
  if (length(n.photos) == 1){
    n.photos <- rep(n.photos, n.animals)
  } else {
    if (length(n.photos) != n.animals){
      stop("The length of 'n.photos' should be equal to 'n.animals'.")
    }
  }

  pars <- c(mus, sigmas, rhos, psis, phis)

  # Organise parameters into useful components
  par.list <- organise.pars(pars, n.animals, n.photos, m, block.only = TRUE)

  # Simulate true measurements for each animal
  true.measurements <- rmvnorm(n.animals, par.list$mus, par.list$sigma)

  # Simulate drone measurements for each photo
  drone.measurements <- vector("list", n.animals)
  for (i in 1:n.animals){
    drone.measurements[[i]] <- rmvnorm(n.photos[i],
                                       true.measurements[i, ],
                                       par.list$xi)
  }

  # Output data frame
  animal.id <- rep(rep(1:n.animals, n.photos), each = m)
  photo.id <- unlist(lapply(lapply(n.photos,
                                   function(x) 1:x),
                            function(x) rep(x, each = m)))
  dim <- rep(1:m, sum(n.photos))

  data.frame(
    animal.id = factor(animal.id),
    photo.id = factor(photo.id),
    dim = factor(dim),
    measurement = unlist(lapply(drone.measurements, function(x) c(t(x))))
  )
}

# -------------------------------------------------------------------------------------------------------

#' Conduct a Simulation Study for Morphometric Models
#'
#' Simulates multiple datasets and fits a model to each one.
#'
#' @param n.sims Integer. The number of data sets to simulate.
#' @param progressbar Logical. If `TRUE` a progress bar will be
#'     displayed.
#' @param n.cores Integer. The number of cores for parallel processing.
#' @inheritParams sim.measurements
#' @inheritParams fit.morph
#'
#' @return An object of class `lme.morph.sim`. The best way to extract
#'     results from this object is to use [`extract.sim.morph()`]. See
#'     the example below.
#'
#' @seealso [`extract.sim.morph()`] to extract results from the object
#'     returned by this function.
#'
#' @examples
#' ## Running a small simulation study.
#' sim.fits <- sim.morph(n.sims = 5,
#'                       n.animals = 10,
#'                       n.photos = 3,
#'                       mus = c(315, 150, 100),
#'                       sigmas = c(25, 15, 10),
#'                       rhos = c(0.85, 0.80, 0.75),
#'                       psis = c(10, 6, 4),
#'                       phis = c(0.5, 0.4, 0.3),
#'                       method = "REML",
#'                       progressbar = FALSE,
#'                       n.cores = 1)
#' ## Extracting p-values from tests for isometry from each model
#' ## fit.
#' iso.p <- function(x) summary(x, type = "isometric-pca")[, 3]
#' extract.sim.morph(sim.fits, FUN = iso.p)
#'
#' @export
sim.morph <- function(n.sims, n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                      method = "REML", progressbar = TRUE, n.cores = 1){
  # If n.photos is scalar, then apply it to all individuals
  if (length(n.photos) == 1){
    n.photos <- rep(n.photos, n.animals)
  } else {
    if (length(n.photos) != n.animals){
      stop("The length of 'n.photos' should be equal to 'n.animals'.")
    }
  }

  # Number of dimensions
  m <- length(mus)

  # Parameter validation
  if (length(sigmas) != m){
    stop("The 'sigmas' argument must have an element for each dimension.")
  }
  if (length(psis) != m){
    stop("The 'psis' argument must have an element for each dimension.")
  }
  if (length(rhos) != sum(1:(m - 1))){
    stop("The 'rhos' argument must have an element for each pair of dimensions.")
  }
  if (length(phis) != sum(1:(m - 1))){
    stop("The 'phis' argument must have an element for each pair of dimensions.")
  }

  # Combine parameters into vector
  pars <- c(mus, sigmas, rhos, psis, phis)

  # Function for parallel processing
  sim_one <- function(x, n.animals, n.photos, m, pars, method){
    data <- sim.measurements(n.animals, n.photos, m, mus, sigmas, rhos, psis, phis)
    try(fit.morph(data, method = method), silent = TRUE)
  }

  # Running simulations
  if (n.cores == 1){
    # Sequential processing with optional progress bar
    fits <- vector(mode = "list", length = n.sims)
    if (progressbar){
      pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
    }
    for (i in 1:n.sims){
      fits[[i]] <- sim_one(i, n.animals, n.photos, m, pars, method)
      if (progressbar){
        setTxtProgressBar(pb, i)
      }
    }
    if (progressbar){
      close(pb)
    }
  } else {
    # Paralllel processing
    cl <- makeCluster(n.cores)
    on.exit(stopCluster(cl))
    fits <- pblapply(1:n.sims, sim_one, n.animals, n.photos, m, pars,
                     method, cl = cl)
  }

  # Preparing output
  settings <- list(
    n.sims = n.sims,
    n.animals = n.animals,
    n.photos = n.photos,
    mus = mus,
    sigmas = sigmas,
    rhos = rhos,
    psis = psis,
    phis = phis,
    method = method
  )
  converged <- sapply(fits, function(x) class(x)[1] != "try-error")
  n.not.converged <- sum(!converged)
  if (n.not.converged > 0){
      warning(paste0("A total of ", n.not.converged, " model ", ifelse(n.not.converged == 1, "fit", "fits"), " did not converge and ", ifelse(n.not.converged == 1, "has", "have"), " been omitted from the simulation results."))
      fits <- fits[converged]
  }
  out <- list(settings = settings, fits = fits)
  class(out) <- "lme.morph.sim"
  out
}

# -------------------------------------------------------------------------------------------------------

#' Extract Results from a Morphometric Simulation Study
#'
#' Extracts results from each model fitted in a simulation study
#' conducted using [`sim.morph()`].
#'
#' @param sim.res An object of class `lme.morph.sim`, returned by
#'     [`sim.morph()`]
#' @param FUN A function to apply to each model fit. Defaults to
#'     [`summary.lme.morph()`], which extracts parameter estimates and
#'     standard errors.
#' @inheritParams sim.morph
#'
#' @return An array containing the results of applying FUN to each
#'     model fit.
#'
#' @inherit sim.morph examples
#'
#' @export
extract.sim.morph <- function(sim.res, FUN = summary, n.cores = 1){
  if (!inherits(sim.res, "lme.morph.sim")) {
    stop("'sim.res' must be an object of class 'lme.morph.sim'")
  }

  if (n.cores > 1) {
    # Parallel processing
    cl <- makeCluster(n.cores)
    on.exit(stopCluster(cl))

    # Export necesssary functions
    clusterExport(cl, "FUN", envir = environment())

    # Apply function to all fits
    out <- simplify2array(pblapply(sim.res$fits, FUN, cl = cl))
  } else {
    # Sequential processing
    out <- simplify2array(lapply(sim.res$fits, FUN))
  }

  out
}
