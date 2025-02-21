# -------------------------------------------------------------------------------------------------------

# calc.conditional.ratio

test_that("calc.conditional.ratio produces expected results", {
  # Test data
  test_data <- sim.measurements(
    n.animals = 50,
    n.photos = rep(10, 50),
    m = 3,
    pars = c(315, 150, 100,    # means
             25, 15, 10,       # SDs
             0.85, 0.80, 0.75,  # correlations
             10, 6, 4,    # measurement SDs
             0.5, 0.4, 0.3)  # measurement correlations
  )

  # Fit model to test data
  fit <- fit.morph(test_data)

  # Test basic functionality for lm
  newdata <- seq(100, 300, length.out = 5)
  result_lm <- calc.conditional.ratio(fit, y.dim = 2, x.dim = 1,
                                      newdata.x.dim = newdata, type = "lm")

  # Check output structure
  expect_true(is.matrix(result_lm))
  expect_equal(colnames(result_lm), c("Estimate", "Std. Error"))
  expect_equal(nrow(result_lm), length(newdata))

  # Check vals are reasonable
  # Ratios should be pos:
  expect_true(all(result_lm[, "Estimate"] > 0))
  # ses should be pos:
  expect_true(all(result_lm[, "Std. Error"] > 0))

  # Test PCA type
  result_pca <- calc.conditional.ratio(fit, y.dim = 2, x.dim = 1,
                                       newdata.x.dim = newdata, type = "pca")

  # Check outut structure
  expect_true(is.matrix(result_pca))
  expect_equal(colnames(result_pca), c("Estimate", "Std. Error"))
  expect_equal(nrow(result_pca), length(newdata))

  # Check vals are reasonable
  expect_true(all(result_pca[, "Estimate"] > 0))
  expect_true(all(result_pca[, "Std. Error"] > 0))

  # lm and PCA should give dif results
  expect_false(identical(result_lm[, "Estimate"], result_pca[, "Estimate"]))

  # Test with dif dims
  result_23 <- calc.conditional.ratio(fit, y.dim = 2, x.dim = 3,
                                      newdata.x.dim = newdata[1:2], type = "lm")
  expect_equal(nrow(result_23), 2)
  expect_true(all(result_23[, "Estimate"] > 0))

  # Test error handling
  expect_error(
    calc.conditional.ratio(fit, y.dim = 1, x.dim = 4, newdata.x.dim = newdata),
    "Dimension index exceeds number of dimensions in data"
  )

  expect_error(
    calc.conditional.ratio(fit, y.dim = 1, x.dim = 1, newdata.x.dim = newdata),
    "Numerator and denominator dimensions must be different"
  )
})

test_that("calc.conditional.ratio handles edge cases", {
  # Test data w/conservative values
  test_data <- sim.measurements(
    n.animals = 5,         # More animals for better convergence
    n.photos = rep(3, 5),  # More pics per animal
    m = 3,
    pars = c(150, 100, 50, #Smallee means
             30, 20, 10,  # Smaller SDs
             0.5, 0.5, 0.5, # Correlations
             1.2, 1.0, 0.8, # Much smaller measurement SDs
             0.2, 0.2, 0.2) # Smaller measurement correlations
  )

    fit <- fit.morph(test_data)

  # Test w/moderate x vals
  test_x <- c(75, 150)
  result <- calc.conditional.ratio(fit, y.dim = 2, x.dim = 1,
                                   newdata.x.dim = test_x, type = "lm")

  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))

  # Test w/single val
  single_x <- 100
  result_single <- calc.conditional.ratio(fit, y.dim = 2, x.dim = 1,
                                          newdata.x.dim = single_x, type = "lm")
  expect_equal(nrow(result_single), 1)
  expect_true(all(is.finite(result_single)))
})

