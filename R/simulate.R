#' Simulation Functions for morphErr
#'
#' This file contains functions for simulating morphometric data:
#' - sim.measurements(): Simulates measurement data from a survey
#' - sim.morph(): Simulates multiple datasets and fits models
#' - extract.sim.morph(): Extracts results from simulation studies
#'
#' @name simulate
#' @keywords internal
NULL

#' Simulate Morphometric Measurements
#'
#' Simulates data from a morphometric survey, including measurement error.
#'
#' @param n.animals Integer, number of animals to simulate
#' @param n.photos Integer vector specifying number of photos per animal.
#'   If scalar, that number of photos is used for all animals.
#' @param m Integer, number of dimensions to measure
#' @param pars Numeric vector of parameters in the order:
#'   \itemize{
#'     \item mus: mean parameters (length m)
#'     \item sigmas: standard deviations (length m)
#'     \item rhos: correlations (length m*(m-1)/2)
#'     \item psis: measurement error SDs (length m)
#'     \item phis: measurement error correlations (length m*(m-1)/2)
#'   }
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item animal.id: Factor indicating which animal
#'     \item photo.id: Factor indicating which photo
#'     \item dim: Factor indicating which dimension
#'     \item measurement: The simulated measurement value
#'   }
#'
#' @examples
#' # Simulate data for 2 animals, 2 photos each, measuring 2 dimensions
#' pars <- c(100, 50, # means for dimensions 1 and 2
#'           10, 5,  # standard deviations
#'           0.7,  # correlation between dimensions
#'           1, 0.5, # measurement error SDs
#'           0.3)   # measurement error correlation
#' data <- sim.measurements(2, rep(2, 2), 2, pars)
#'
#' @export
sim.measurements <- function(n.animals, n.photos, m, pars){
  # If n.photos is scalar, repliicate for all animals
  if (length(n.photos) == 1){
    n.photos <- rep(n.photos, n.animals)
  } else {
    if (length(n.photos) != n.animals){
      stop("The length of 'n.photos' should be equal to 'n.animals'.")
    }
  }

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

#' Simulation Study for Morphometric Models
#'
#' Simulates multiple datasets and fits models to each one.
#'
#' @param n.sims Integer, number of datasets to simulate
#' @param n.animals Integer, number of animals to simulate
#' @param n.photos Integer vector or scalar, number of photos per animal
#' @param mus Numeric vector of mean parameters for each dimension
#' @param sigmas Numeric vector of standard deviations for each dimension
#' @param rhos Numeric vector of correlations between dimensions
#' @param psis Numeric vector of measurement error SDs
#' @param phis Numeric vector of measurement error correlations
#' @param method Character, either "ML" or "REML"
#' @param progressbar Logical, whether to show a progress bar
#' @param n.cores Integer, number of cores for parallel processing
#'
#' @return An object of class "lme.morph.sim" containing:
#'   \itemize{
#'     \item settings: List of simulation parameters
#'     \item fits: List of fitted models
#'   }
#'
#' @examples
#' \dontrun{
#' # Small simulation study
#' results <- sim.morph(
#'   n.sims = 5,
#'   n.animals = 10,
#'   n.photos = 3,
#'   mus = c(100, 50),
#'   sigmas = c(10, 5),
#'   rhos = 0.7,
#'   psis = c(1, 0.5),
#'   phis = 0.3,
#'   method = "REML",
#'   n.cores = 1
#' )
#' }
#'
#' @export
sim.morph <- function(n.sims, n.animals, n.photos, mus, sigmas, rhos, psis, phis,
                      method = "REML", progressbar = TRUE, n.cores = 1){
  # If n.photos is scalar, then aply it to all individuals
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
    data <- sim.measurements(n.animals, n.photos, m, pars)
    suppressWarnings(fit.morph(data, method = method))
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

  out <- list(settings = settings, fits = fits)
  class(out) <- "lme.morph.sim"
  out
}

#' Extract Results from Simulation Study
#'
#' Extracts specific results from each model fit in a simulation study.
#'
#' @param sim.res An object of class "lme.morph.sim" returned by sim.morph()
#' @param FUN Function to apply to each model fit. Default is summary.
#' @param n.cores Integer, number of cores for parallel processing.
#'
#' @return An array containing the results of applying FUN to each model fit.
#'
#' @examples
#' \dontrun{
#' # Run simulation study
#' results <- sim.morph(
#'   n.sims = 5,
#'   n.animals = 10,
#'   n.photos = 3,
#'   mus = c(100, 50),
#'   sigmas = c(10, 5),
#'   rhos = 0.7,
#'   psis = c(1, 0.5),
#'   phis = 0.3
#' )
#'
#' # Extract summaries
#' summaries <- extract.sim.morph(results)
#' }
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
