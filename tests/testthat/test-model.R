# construct.varcov() test

test_that("construct.varcov works correctly", {
  # Test 1: simple 3x3 block
  sds <- c(25, 15, 10)
  cors <- c(0.85, 0.80, 0.75)
  result <- construct.varcov(sds, cors, 1, 3, block.only = TRUE)
  expect_equal(dim(result), c(3, 3))
  expect_equal(diag(result), sds^2)

  # Test 2: multiple blocks
  expect_equal(
    dim(construct.varcov(sds, cors, 2, 3, block.only = FALSE)),
    c(6, 6)
  )
})

# -------------------------------------------------------------------------------------------------------

# organise.pars() test

test_that("organise.pars works correctly", {
  # Simple test
  pars <- c(315, 150, 100,    # mus
            25, 15, 10,       # sigmas
            0.85, 0.80, 0.75, # rhos
            10, 6, 4,         # psis
            0.5, 0.4, 0.3)    # phis
  result <- organise.pars(pars, 2, c(1,1), 3, block.only = TRUE)

  expect_named(result, c("mus", "sigma", "xi"))
  expect_equal(result$mus, c(315, 150, 100))
})

# -------------------------------------------------------------------------------------------------------

# fit.morph() test

test_that("fit.morph handles inputs correctly", {
  # Reproducibility
  set.seed(123)

  # Measurements
  n_animals <- 10
  n_photos <- 5

  # true measurements for each animal
  true_means <- c(315, 150, 100)
  true_sds <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)

  # Construct correlation matrix
  cor_mat <- matrix(1, 3, 3)
  cor_mat[upper.tri(cor_mat)] <- rhos
  cor_mat[lower.tri(cor_mat)] <- rhos[3:1]

  # Construct covariance matrix
  sigma <- diag(true_sds) %*% cor_mat %*% diag(true_sds)
  true_measurements <- MASS::mvrnorm(n_animals, true_means, sigma)

  # Measurement error for each pic
  measurements <- matrix(0, nrow = n_animals * n_photos, ncol = 3)
  for(i in 1:n_animals) {
    row_idx <- ((i-1)*n_photos + 1):(i*n_photos)
    measurements[row_idx, ] <- matrix(
      rnorm(n_photos * 3,
            rep(true_measurements[i,], each = n_photos),
            # Measurement error SDs
            c(5, 3, 2)),
      ncol = 3
    )
  }

  test_data <- data.frame(
    animal.id = factor(rep(1:n_animals, each = n_photos)),
    photo.id = factor(rep(1:n_photos, n_animals)),
    dim = factor(rep(1:3, each = n_animals * n_photos)),
    measurement = c(measurements)
  )

  # Test basic model fit
  fit <- fit.morph(test_data)

  # Check class
  expect_s3_class(fit, "lme.morph")

  # Check components
  expect_true(is.numeric(fit$coefficients$fixed))
  expect_false(fit$intercept)

  # Test with intercept
  fit_int <- fit.morph(test_data, intercept = TRUE)
  expect_true(fit_int$intercept)
})
