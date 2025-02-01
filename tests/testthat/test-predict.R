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

  # Test invalid object class
  invalid_obj <- list()
  class(invalid_obj) <- "invalid"
  expect_error(
    predict.lme.morph(invalid_obj, y.dim = 1),
    "Invalid model object"
  )

  # Test invalid newdata type
  expect_error(
    predict.lme.morph(fit, y.dim = 1, newdata = matrix(1, 1, 1)),
    "'newdata' must be a data frame"
  )

  # Test invalid type
  expect_error(
    predict.lme.morph(fit, y.dim = 1, newdata = data.frame(dim2 = 130), type = "invalid"),
    "Invalid type argument"
  )
})

test_that("predict.lme.morph handles predictions correctly", {
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

  # Test pred with newdata
  newdata <- data.frame(dim2 = 130)
  result <- predict.lme.morph(fit, y.dim = 1, newdata = newdata, type = "lm")

  # Check output structure
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(nrow(result), 1)

  # Test both pred types
  lm.pred <- predict.lme.morph(fit, y.dim = 1, newdata = newdata, type = "lm")
  pca.pred <- predict.lme.morph(fit, y.dim = 1, newdata = newdata, type = "pca")
  expect_false(identical(lm.pred[1], pca.pred[1]))
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

  # Test pred with no measurements
  result <- predict.lme.morph(fit, y.dim = 1)

  # Shuld return pop means
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(nrow(result), 1)
  expect_true(!is.na(result[1,1]))
  expect_true(abs(result[1,1] - mus[1]) < 50)
})

test_that("predict.lme.morph calculates standard errors correctly", {
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

  # Test w/newdata
  newdata <- data.frame(dim2 = 130)
  result <- predict.lme.morph(fit, y.dim = 1, newdata = newdata)

  # Check se reasoanlble
  expect_false(is.na(result[1,2]))
  expect_true(result[1,2] > 0)
  expect_true(result[1,2] < max(pars[1:3]))
})

# -------------------------------------------------------------------------------------------------------

# predict.from.obs() tests

test_that("predict.from.obs handles input validation correctly", {
  # Setup
  set.seed(1234)
  pars <- c(290, 125, 75, 45, 25, 15, 0.75, 0.80, 0.85, 2.00, 1.50, 1.00, 0.40, 0.50, 0.60)
  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test invalid model object
  invalid_obj <- list()
  class(invalid_obj) <- "invalid"
  expect_error(
    predict.from.obs(invalid_obj, true = c(NA, 130, 75)),
    "Invalid model object"
  )

  # Test invalid true vecter length
  expect_error(
    predict.from.obs(fit, true = c(NA, 130)),
    "Length of 'true' must match"
  )

  # Test no dims to pred
  expect_error(
    predict.from.obs(fit, true = c(290, 130, 75)),
    "No dimensions to predict"
  )
})

test_that("predict.from.obs handles predictions correctly", {
  # Setup
  set.seed(1234)
  pars <- c(290, 125, 75, 45, 25, 15, 0.75, 0.80, 0.85, 2.00, 1.50, 1.00, 0.40, 0.50, 0.60)
  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test basic pred
  true <- c(NA, 130, 75)
  obs <- matrix(c(285, 128, 73), nrow = 1)
  result <- predict.from.obs(fit, true = true, obs = obs)

  # Check output structure
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)
  # Pred val not NA:
  expect_false(is.na(result[1]))
  # Known values unchanged:
  expect_equal(result[2:3], true[2:3])
})
