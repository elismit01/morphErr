# construct.varcov() test

test_that("construct.varcov works correctly", {
  # Test 1: simple 2x2 block
  sds <- c(1, 2)
  cors <- 0.5
  expect_equal(
    construct.varcov(sds, cors, 1, 2, block.only = TRUE),
    matrix(c(1, 1, 1, 4), 2, 2)
  )

  # Test 2: multiple blocks
  expect_equal(
    dim(construct.varcov(sds, cors, 2, 2, block.only = FALSE)),
    c(4, 4)
  )
})

# -------------------------------------------------------------------------------------------------------

# organise.pars() test

test_that("organise.pars works correctly", {
  # Simple test
  pars <- c(10, 20,    # mus
            2, 3,      # sigmas
            0.5,        # rhos
            0.1, 0.2,  # psis
            0.3)       # phis
  result <- organise.pars(pars, 2, c(1,1), 2, block.only = TRUE)

  expect_named(result, c("mus", "sigma", "xi"))
  expect_equal(result$mus, c(10, 20))
})

# -------------------------------------------------------------------------------------------------------

# fit.morph() test

test_that("fit.morph handles inputs correctly", {
  # Reproducibility
  set.seed(123)

  # Measurements
  n_animals <- 5
  n_photos <- 3
  rho <- 0.7

  # true measurements for each animal
  true_means <- c(100, 50)
  true_sds <- c(10, 5)
  sigma <- matrix(
    c(true_sds[1]^2, rho*prod(true_sds),
      rho*prod(true_sds), true_sds[2]^2),
    nrow = 2
  )
  true_measurements <- MASS::mvrnorm(n_animals, true_means, sigma)

  # Measurement error for each pic
  measurements <- matrix(0, nrow = n_animals * n_photos, ncol = 2)
  for(i in 1:n_animals) {
    row_idx <- ((i-1)*n_photos + 1):(i*n_photos)
    measurements[row_idx, ] <- matrix(
      rnorm(n_photos * 2,
            rep(true_measurements[i,], each = n_photos),
            # Measurement error SDs
            c(1, 0.5)),
      ncol = 2
    )
  }

  test_data <- data.frame(
    animal.id = factor(rep(1:n_animals, each = n_photos)),
    photo.id = factor(rep(1:n_photos, n_animals)),
    dim = factor(rep(1:2, each = n_animals * n_photos)),
    measurement = c(measurements)
  )

  # Test basic model fit
  fit <- suppressWarnings(fit.morph(test_data))

  # Check class
  expect_s3_class(fit, "lme.morph")

  # Check components
  expect_true(is.numeric(fit$coefficients$fixed))
  expect_false(fit$intercept)

  # Test with intercept
  fit_int <- suppressWarnings(fit.morph(test_data, intercept = TRUE))
  expect_true(fit_int$intercept)
})
