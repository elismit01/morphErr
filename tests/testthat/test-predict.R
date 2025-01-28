# calc.betas() tests

test_that("calc.betas handles input validation", {
  # Test data
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test PCA with multiple predicters
  expect_error(
    calc.betas(fit, y.dim = 1, x.dim = c(2, 3), type = "pca"),
    "Type 'pca' is only available for single predictor relationships"
  )

  # Test output structure
  result <- calc.betas(fit, y.dim = 1, x.dim = 2, type = "lm")
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(rownames(result), c("beta0", "beta2"))
})

test_that("calc.betas handles different types correctly", {
  # Test data setup (using the same)
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test lm type
  lm_result <- calc.betas(fit, y.dim = 1, x.dim = 2, type = "lm")
  expect_true(!is.null(lm_result))
  expect_equal(nrow(lm_result), 2)

  # Test PCa type
  pca_result <- calc.betas(fit, y.dim = 1, x.dim = 2, type = "pca")
  expect_true(!is.null(pca_result))
  expect_equal(nrow(pca_result), 2)
})

# -------------------------------------------------------------------------------------------------------

# predict.lme.morph() tests

test_that("predict.lme.morph handles input validation correctly", {
  # Data set up
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test invalid object class
  invalid_obj <- list()
  class(invalid_obj) <- "invalid"
  expect_error(
    predict.lme.morph(invalid_obj, y.dim = 1),
    "Object must be of class 'lme.morph'"
  )

  # Test providing both measurement types
  expect_error(
    predict.lme.morph(fit,
                      true.measurements = data.frame(dim1 = 1),
                      observed.measurements = matrix(1, 1, 1),
                      y.dim = 1),
    "Cannot provide both true and observed measurements"
  )

  # Test invalid true.measurements type
  expect_error(
    predict.lme.morph(fit,
                      true.measurements = matrix(1, 1, 1),
                      y.dim = 1),
    "true.measurements must be a data frame"
  )
})

test_that("predict.lme.morph handles true measurements correctly", {
  # Data setup
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test prediction wit true measurements
  true.meas <- data.frame(dim2 = 130)
  result <- predict.lme.morph(fit,
                              true.measurements = true.meas,
                              y.dim = 1,
                              type = "lm")

  # Check outpt structure
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(nrow(result), 1)

  # Test both pred types
  lm.pred <- predict.lme.morph(fit,
                               true.measurements = true.meas,
                               y.dim = 1,
                               type = "lm")
  pca.pred <- predict.lme.morph(fit,
                                true.measurements = true.meas,
                                y.dim = 1,
                                type = "pca")

  expect_false(identical(lm.pred[1], pca.pred[1]))
})

test_that("predict.lme.morph handles observed measurements correctly", {
  # Setup again
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test with a single relistic obs
  obs.mat <- matrix(c(285, 128, 73), nrow = 1)
  result.mat <- suppressWarnings(
    predict.lme.morph(fit,
                      observed.measurements = obs.mat,
                      y.dim = 1)
  )

  # Check basic structure + values
  expect_true(is.matrix(result.mat))
  expect_equal(colnames(result.mat), c("Estimate", "Std. Error"))
  expect_true(!is.na(result.mat[1,1]))
  expect_true(is.na(result.mat[1,2]))

  # Test that prediction is somewhat close to the observation
  expect_true(abs(result.mat[1,1] - obs.mat[1,1]) < 50)
})

test_that("predict.lme.morph handles no measurements case correctly", {
  # Setup
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test prediction with no measurements
  result <- predict.lme.morph(fit, y.dim = 1)

  # Should return pop means
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(nrow(result), 1)
  expect_true(!is.na(result[1,1]))
  # (Noot testing SE anymore cause it might be NA)
  expect_true(abs(result[1,1] - mus[1]) < 50)
})

test_that("predict.lme.morph calculates standard errors for true measurements", {
  # Setup
  set.seed(1234)
  mus <- c(290, 125, 75)
  sigmas <- c(45, 25, 15)
  rhos <- c(0.75, 0.80, 0.85)
  psis <- c(2.00, 1.50, 1.00)
  phis <- c(0.40, 0.50, 0.60)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test w/true measurements
  true.meas <- data.frame(dim2 = 130)
  result <- predict.lme.morph(fit, true.measurements = true.meas, y.dim = 1)

  # Check se is calculated and reasonable
  expect_false(is.na(result[1,2]))
  expect_true(result[1,2] > 0)
  # Compar to max of mean dimensions (cause se shouldn't be larger than the measurement scale):
  expect_true(result[1,2] < max(pars[1:3]))
})
