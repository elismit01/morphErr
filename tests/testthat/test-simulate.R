# sim.measurements() tests

test_that("sim.measurements handles inputs correctly", {
  # Test 1: basic functionality with equal photos
  pars <- c(100, 50,  # means
            10, 5,  # SDs
            0.7,  # correlation
            1, 0.5,  # measurement SDs
            0.3) # measurement correlation

  # 3 animals, 2 photos each:
  data1 <- sim.measurements(3, 2, 2, pars)

  # Basic structure checkss
  expect_s3_class(data1, "data.frame")
  expect_named(data1, c("animal.id", "photo.id", "dim", "measurement"))
  expect_equal(nrow(data1), 3 * 2 * 2)  # n_animals * n_photos * n_dims

  # Test 2: dif numbers of photos per animal
  n_photos <- c(2, 3, 1)
  data2 <- sim.measurements(3, n_photos, 2, pars)
  expect_equal(nrow(data2), sum(n_photos) * 2)  # total_photos * n_dims

  # Test 3: error handling
  expect_error(
    sim.measurements(2, c(2, 2, 2), 2, pars),
    "length of 'n.photos' should be equal to 'n.animals'"
  )
})

# -------------------------------------------------------------------------------------------------------

# sim.morph() tests

test_that("sim.morph runs simulations correctly", {
  # Set up small simulation
  n.sims <- 2
  n.animals <- 3
  n.photos <- 2
  mus <- c(100, 50)
  sigmas <- c(10, 5)
  rhos <- 0.7
  psis <- c(1, 0.5)
  phis <- 0.3

  # Run simulation
  set.seed(123)
  result <- sim.morph(
    n.sims = n.sims,
    n.animals = n.animals,
    n.photos = n.photos,
    mus = mus,
    sigmas = sigmas,
    rhos = rhos,
    psis = psis,
    phis = phis,
    progressbar = FALSE
  )

  # check structure
  expect_s3_class(result, "lme.morph.sim")
  expect_named(result, c("settings", "fits"))
  expect_equal(length(result$fits), n.sims)

  # Check settings are stored correctly
  expect_equal(result$settings$n.sims, n.sims)
  expect_equal(result$settings$mus, mus)

  # Check model fits
  expect_s3_class(result$fits[[1]], "lme.morph")

  # Test error handling
  expect_error(
    sim.morph(n.sims = 2, n.animals = 3, n.photos = 2,
              mus = c(100, 50), sigmas = c(10), rhos = 0.7,
              psis = c(1, 0.5), phis = 0.3),
    "The 'sigmas' argument must have an element for each dimension"
  )
})

# -------------------------------------------------------------------------------------------------------

# extract.sim.morph() tests

test_that("extract.sim.morph works correctly", {
  # Create a small simulation study
  set.seed(123)
  sim_result <- sim.morph(
    n.sims = 2,
    n.animals = 3,
    n.photos = 2,
    mus = c(100, 50),
    sigmas = c(10, 5),
    rhos = 0.7,
    psis = c(1, 0.5),
    phis = 0.3,
    progressbar = FALSE
  )

  # Test basic extraction
  extracted <- extract.sim.morph(sim_result)
  expect_true(is.array(extracted))
  # Should be at least 2D
  expect_true(length(dim(extracted)) >= 2)

  # Test custom function
  get_fixed <- function(fit) fit$coefficients$fixed
  fixed_effects <- extract.sim.morph(sim_result, FUN = get_fixed)
  expect_true(is.array(fixed_effects))
  # Number of simulations
  expect_equal(dim(fixed_effects)[2], 2)

  # Test eror handling
  expect_error(
    extract.sim.morph(list(a = 1)),
    "'sim.res' must be an object of class 'lme.morph.sim'"
  )
})
