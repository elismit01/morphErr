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
#'     2, and so on, where `eqn{m}` is the number of dimensions.
#'
#' For example, if \eqn{m = 4}, then the first three elements are the
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
#' @param data A data frame with columns `animal.id`, `photo.id`, and
#'     `dim`, provided instead of `n.animals` and `n.photos`. This
#'     provides the user with full control over which dimensions are
#'     measured in which photos from which animals.
#' @param mus A vector with an element for each dimension, providing
#'     the means of the true dimension sizes in the population.
#' @param sigmas A vector with an element for each dimension,
#'     providing the standard deviations for true dimension sizes in
#'     the population.
#' @param rhos A vector, with one element for each pair of dimensions,
#'     providing the pairwise correlations between true dimension
#'     sizes in the population. See 'Details' for the correct order
#'     for the correlations.
#' @param psis A vector with an element for each dimension, providing
#'     the standard deviations of measurement errors for the
#'     dimensions.
#' @param phis A vector, with one element for each pair of dimensions,
#'     providing the pairwise correlations between measurement errors
#'     for the dimensions. See 'Details' for the correct order for the
#'     correlations.
#' @param log.transform Logical. If `TRUE`, the parameters are
#'     considered to correspond to a model where the response was
#'     log-transformed. The data frame returned by this function will
#'     contain the back-transformed measurments.
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
#' sim.data <- sim.measurements(n.animals = 10, n.photos = 2,
#'                              mus = c(315, 150, 100),
#'                              sigmas = c(25, 15, 10),
#'                              rhos = c(0.85, 0.80, 0.75),
#'                              psis = c(10, 6, 4),
#'                              phis = c(0.5, 0.4, 0.3))
#' head(sim.data)
#' 
#' ## Simulating data for two animals, with different numbers of
#' ## photos for each, and different measurements available from
#' ## different photos.
#' data <- data.frame(animal.id = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
#'                    photo.id = c(1, 1, 2, 2, 2, 1, 2, 3, 3, 4, 4),
#'                    dim = c(1, 3, 1, 2, 3, 1, 1, 1, 3, 2, 3))
#' sim.data <- sim.measurements(data = data,
#'                              mus = c(315, 150, 100),
#'                              sigmas = c(25, 15, 10),
#'                              rhos = c(0.85, 0.80, 0.75),
#'                              psis = c(10, 6, 4),
#'                              phis = c(0.5, 0.4, 0.3))
#' sim.data
#'
#' @export
sim.measurements <- function(n.animals = NULL, n.photos = NULL,
                             data = NULL, mus, sigmas, rhos, psis,
                             phis, log.transform = FALSE){
    ## Number of dimensions.
    m <- length(mus)
    ## Sorting out data structure.
    if (is.null(data)){
        if (is.null(n.animals) | is.null(n.photos)){
            stop("Provide either 'data', or both 'n.animals' and 'n.photos'")
        } else {
            ## If n.photos is scalar, replicate for all animals
            if (length(n.photos) == 1){
                n.photos <- rep(n.photos, n.animals)
            } else {
                if (length(n.photos) != n.animals){
                    stop("The length of 'n.photos' should be equal to 'n.animals'.")
                }
            }
            data <- data.frame(animal.id = rep(rep(1:n.animals, times = n.photos), each = m),
                               photo.id = rep(sequence(n.photos), each = m),
                               dim = rep(1:m, sum(n.photos)))
        }
    } else {
        ## Overwriting animal IDs to start at 1 and increment.
        data$animal.id <- match(data$animal.id, unique(data$animal.id))
        ## Total number of animals.
        n.animals <- length(unique(data$animal.id))
        ## Overwriting photo IDs to start at 1 and increment for each animal.
        for (i in 1:n.animals){
            data$photo.id[data$animal.id == i] <- match(data$photo.id[data$animal.id == i],
                                                        unique(data$photo.id[data$animal.id == i]))
        }
        ## Total number of photographs per animal.
        n.photos <- aggregate(photo.id ~ animal.id, data = data, FUN = function(x) length(unique(x)))$photo.id
    }
    pars <- c(mus, sigmas, rhos, psis, phis)
    
    ## Organise parameters into useful components
    par.list <- organise.pars(pars, n.animals, n.photos, m, block.only = TRUE)
    
    ## Simulate true measurements for each animal
    true.measurements <- rmvnorm(n.animals, par.list$mus, par.list$sigma)
    
    ## Simulate drone measurements for each photo
    drone.measurements <- vector("list", n.animals)
    for (i in 1:n.animals){
        drone.measurements[[i]] <- rmvnorm(n.photos[i],
                                           true.measurements[i, ],
                                           par.list$xi)
    }
    
    ## Creating full data frame.
    full.animal.id <- rep(rep(1:n.animals, n.photos), each = m)
    full.photo.id <- unlist(lapply(lapply(n.photos,
                                          function(x) 1:x),
                                   function(x) rep(x, each = m)))
    full.dim <- rep(1:m, sum(n.photos))
    full.measurement <- unlist(lapply(drone.measurements, function(x) c(t(x))))
    full.df <- data.frame(animal.id = full.animal.id, photo.id = full.photo.id,
                          dim = full.dim, measurement = full.measurement)
    ## Checking which rows are required.
    keep <- interaction(full.df$animal.id, full.df$photo.id, full.df$dim) %in%
        interaction(data$animal.id, data$photo.id, data$dim)
    ## Log-transforming measurements if required.
    if (log.transform){
        full.measurement <- exp(full.measurement)
    }
    ## Putting together the full data frame.
    out <- cbind(data, measurement = full.measurement[keep])
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
                      log.transform = FALSE, method = "REML", progressbar = TRUE,
                      n.cores = 1){
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
  sim_one <- function(x, n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                      log.transform, method){
    data <- sim.measurements(n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                             log.transform = log.transform)
    try(fit.morph(data, log.transform = log.transform, method = method), silent = TRUE)
  }

  # Running simulations
  if (n.cores == 1){
    # Sequential processing with optional progress bar
    fits <- vector(mode = "list", length = n.sims)
    if (progressbar){
      pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
    }
    for (i in 1:n.sims){
      fits[[i]] <- sim_one(i, n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                           log.transform, method)
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
    fits <- pblapply(1:n.sims, sim_one, n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                     log.transform, method, cl = cl)
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
    log.transform = log.transform,
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
