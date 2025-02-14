# calc.betas() tests

test_that("calc.betas handles input validation", {
  # Test data
  set.seed(1234)
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
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

test_that("predict.lme.morph handles multiple predictions correctly", {
  # Setup
  set.seed(1234)
  mus <- c(315, 150, 100)
  sigmas <- c(25, 15, 10)
  rhos <- c(0.85, 0.80, 0.75)
  psis <- c(10, 6, 4)
  phis <- c(0.5, 0.4, 0.3)
  pars <- c(mus, sigmas, rhos, psis, phis)

  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Create multiple row newdata (6 vals each)
  newdata <- data.frame(
    dim2 = seq(120, 130, by = 2),
    dim3 = seq(70, 80, by = 2)
  )

  # Test preds
  result <- predict.lme.morph(fit, y.dim = 1, newdata = newdata)

  # Check structure
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("Estimate", "Std. Error"))
  expect_equal(nrow(result), nrow(newdata))

  # Check preds increase w/increasing predictors
  expect_true(all(diff(result[,"Estimate"]) > 0))

  # Check ses are reasonable
  expect_true(all(result[,"Std. Error"] > 0))
  expect_true(all(result[,"Std. Error"] < max(mus) * 2))
})

test_that("predict.lme.morph handles edge cases in multiple predictions", {
  # Setup
  set.seed(1234)
  pars = c(315, 150, 100,    # means
           25, 15, 10,       # SDs
           0.85, 0.80, 0.75,  # correlations
           10, 6, 4,    # measurement SDs
           0.5, 0.4, 0.3)  # measurement correlations
  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test single row newdata works same
  single_newdata <- data.frame(dim2 = 130)
  single_result <- predict.lme.morph(fit, y.dim = 1, newdata = single_newdata)
  expect_equal(nrow(single_result), 1)

  # Test large newdata
  large_newdata <- data.frame(
    dim2 = seq(100, 200, length.out = 100),
    dim3 = seq(50, 100, length.out = 100)
  )
  large_result <- predict.lme.morph(fit, y.dim = 1, newdata = large_newdata)
  expect_equal(nrow(large_result), 100)

  # Test with na values (should error)
  na_data <- data.frame(dim2 = c(130, NA, 140))
  expect_error(
    predict.lme.morph(fit, y.dim = 1, newdata = na_data),
    "missing values"
  )
})

# -------------------------------------------------------------------------------------------------------

# predict.from.obs() tests

test_that("predict.from.obs handles input validation correctly", {
  # Setup
  set.seed(1234)
  pars = c(315, 150, 100,    # means
           25, 15, 10,       # SDs
           0.85, 0.80, 0.75,  # correlations
           10, 6, 4,    # measurement SDs
           0.5, 0.4, 0.3)  # measurement correlations
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
    predict.from.obs(fit, true = c(315, 150, 100)),
    "No dimensions to predict"
  )
})

test_that("predict.from.obs handles predictions correctly", {
  # Setup
  set.seed(1234)
  pars = c(315, 150, 100,    # means
           25, 15, 10,       # SDs
           0.85, 0.80, 0.75,  # correlations
           10, 6, 4,    # measurement SDs
           0.5, 0.4, 0.3)  # measurement correlations
  data <- sim.measurements(10, rep(5, 10), 3, pars)
  fit <- fit.morph(data)

  # Test with named params
  true <- c(NA, 130, 75)
  obs <- matrix(c(285, 128, 73), nrow = 1)
  result1 <- predict.from.obs(fit, true = true, obs = obs)

  # Test with unnamed params
  result2 <- predict.from.obs(fit, true, obs)

  # Check both method give same result
  expect_equal(result1, result2)

  # Check output structure
  expect_true(is.numeric(result1))
  expect_equal(length(result1), 3)
  # Pred val not NA:
  expect_false(is.na(result1[1]))
  # Known values unchanged:
  expect_equal(result1[2:3], true[2:3])
})
