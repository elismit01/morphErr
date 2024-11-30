# mean.rcnorm.approx()

test_that("mean.rcnorm.approx produces expected results", {
  # Basic case with no correlation (allowing for approximation error)
  expect_equal(
    mean.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0),
    1,
    tolerance = 0.02
  )

  # Test with correlation
  result <- mean.rcnorm.approx(mean1 = 2, mean2 = 1, sd1 = 0.2, sd2 = 0.1, rho = 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(!is.na(result))

  # test scaling property (allowing for approximation error)
  ratio1 <- mean.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  ratio2 <- mean.rcnorm.approx(mean1 = 2, mean2 = 2, sd1 = 0.2, sd2 = 0.2, rho = 0)
  expect_equal(ratio1, ratio2, tolerance = 0.02)

  # Test expectd relationship between inputs
  # When mean1 = 2*mean2, expect ratio around 2
  result_double <- mean.rcnorm.approx(mean1 = 2, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  expect_equal(result_double, 2, tolerance = 0.02)
})

# sd.rcnorm.approx()

test_that("sd.rcnorm.approx produces expected results", {
  # First test group: basic output properties
  result <- sd.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result > 0)

  # Second test group: scaling behaviour
  result1 <- sd.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  result2 <- sd.rcnorm.approx(mean1 = 2, mean2 = 2, sd1 = 0.2, sd2 = 0.2, rho = 0)
  expect_equal(result1, result2, tolerance = 0.02)

  # Third test group: corelation effects
  result_uncorr <- sd.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  result_corr <- sd.rcnorm.approx(mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0.5)
  expect_type(result_corr, "double")
  expect_true(result_corr > 0)
})

# drcnorm.approx()

test_that("drcnorm.approx produces expected results", {
  # First test group: basic density behaviour
  w <- seq(0.5, 1.5, by = 0.1)
  result <- drcnorm.approx(w, mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  expect_type(result, "double")
  expect_length(result, length(w))
  expect_true(all(result >= 0))

  # Second test group: expected behaviours with standardised inputs
  w_test <- c(-0.1, 0, 0.1)
  result_std <- drcnorm.approx(w_test, mean1 = 0, mean2 = 1, sd1 = 1, sd2 = 1, rho = 0)
  expect_true(all(result_std >= 0))
  expect_length(result_std, length(w_test))

  # Third test group: corelation behaviour
  result_uncorr <- drcnorm.approx(1, mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  result_corr <- drcnorm.approx(1, mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0.5)
  expect_false(identical(result_uncorr, result_corr))
})

# drcnorm()

test_that("drcnorm validates inputs correctly", {
  expect_error(
    drcnorm("x", mean1 = 0, mean2 = 1, sd1 = 1, sd2 = 1, rho = 0),
    "All arguments must be numeric"
  )

  expect_error(
    drcnorm(1, mean1 = 0, mean2 = 1, sd1 = -1, sd2 = 1, rho = 0),
    "Standard deviations must be positive"
  )

  expect_error(
    drcnorm(1, mean1 = 0, mean2 = 1, sd1 = 1, sd2 = 1, rho = 1.5),
    "Correlation coefficient must be between -1 and 1"
  )
})

test_that("drcnorm produces expected values", {
  # Test output type and length
  w <- seq(0.5, 1.5, by = 0.1)
  result <- drcnorm(w, mean1 = 0, mean2 = 1, sd1 = 1, sd2 = 1, rho = 0)
  expect_type(result, "double")
  expect_length(result, length(w))

  # Test all outputs are non-negative
  expect_true(all(result >= 0))

  # Test symetry for standard case
  w_sym <- c(-1, 1)
  result_sym <- drcnorm(w_sym, mean1 = 0, mean2 = 1, sd1 = 1, sd2 = 1, rho = 0)
  expect_equal(result_sym[1], result_sym[2])

  # Test consistency with drcnorm.approx for simple cases
  w_test <- 1
  exact <- drcnorm(w_test, mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  approx <- drcnorm.approx(w_test, mean1 = 1, mean2 = 1, sd1 = 0.1, sd2 = 0.1, rho = 0)
  expect_equal(exact, approx, tolerance = 0.1)
})

