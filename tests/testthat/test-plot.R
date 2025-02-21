# Create data once for all tests
create_test_data <- function() {
  set.seed(123)
  sim.measurements(
    n.animals = 10,
    n.photos = 3,
    m = 3,
    pars = c(315, 150, 100,    # means
             25, 15, 10,       # SDs
             0.85, 0.80, 0.75,  # correlations
             10, 6, 4,    # measurement SDs
             0.5, 0.4, 0.3)  # measurement correlations
  )
}

# -------------------------------------------------------------------------------------------------------

# plot.morph()

test_that("plot.morph creates correct plot structure", {
  # Get test data
  test_data <- create_test_data()

  # Test basic plotting
  expect_silent(plot.morph(test_data))

  # test ratio plotting
  expect_silent(plot.morph(test_data, ratios = TRUE))

  # Test custom limits and labels
  expect_silent(plot.morph(test_data,
                           xlim = c(0, 200),
                           ylim = c(0, 100),
                           xlab = "Width",
                           ylab = "Height"))

  # Test plot.data = FALSE
  expect_silent(plot.morph(test_data, plot.data = FALSE))
})

test_that("plot.morph handles invalid inputs", {
  # Create minimal test data with only one dimension
  bad_data <- data.frame(
    animal.id = factor(1:2),
    photo.id = factor(1:2),
    dim = factor(1),
    measurement = 1:2
  )

  # Should errror when requesting non-existent dimension
  expect_error(
    plot.morph(bad_data, dims = c(1, 2)),
    "Not all requested dimensions are present in the data"
  )
})

test_that("plot.morph handles edge cases", {
  test_data <- create_test_data()

  # Test with single animal
  single_animal <- test_data[test_data$animal.id == levels(test_data$animal.id)[1], ]
  expect_silent(plot.morph(single_animal))

  # Test with single photo per animal
  single_photo <- test_data[test_data$photo.id == levels(test_data$photo.id)[1], ]
  expect_silent(plot.morph(single_photo))
})

test_that("plot.lme.morph creates basic plots correctly", {
  # Create and fit test data
  test_data <- create_test_data()
  fit <- fit.morph(test_data)

  # Test different plot types
  expect_silent(plot(fit, type = "data"))
  expect_silent(plot(fit, type = "ratio"))
})

# -------------------------------------------------------------------------------------------------------

# plot.lme.morph()

test_that("plot.lme.morph handles line overlays", {
  # Create and fit test data
  test_data <- create_test_data()
  fit <- fit.morph(test_data)

  # Test different line types
  expect_silent(plot(fit, type = "data", line.type = "lm"))
  expect_silent(plot(fit, type = "data", line.type = "pca"))

  # Test with cis
  expect_silent(plot(fit, type = "data", line.type = "lm", confints = TRUE))
})

test_that("plot.lme.morph handles adding to existing plots", {
  # Create + fit test data
  test_data <- create_test_data()
  fit <- fit.morph(test_data)

  # Test first plot
  expect_silent(plot(fit, type = "data", line.type = "lm"))

  # Test adding to existing plot
  expect_silent(plot(fit, type = "data", line.type = "pca", add = TRUE, lty = 2))
})

test_that("plot.lme.morph handles different dimension combinations correctly", {
  # Test data
  test_data <- create_test_data()
  fit <- fit.morph(test_data)

  # Test diff dim combinations
  for (type in c("pca", "lm")) {
    # Test plots
    expect_silent(plot(fit, dims = c(1, 2), line.type = type))
    expect_silent(plot(fit, dims = c(1, 3), line.type = type))
    expect_silent(plot(fit, dims = c(2, 3), line.type = type))

    # Get preds for each combination to varify dif dims are used
    x_val <- 200  # (test value)
    pred_12 <- predict(fit, y.dim = 2, newdata = data.frame(dim1 = x_val), type = type)
    pred_13 <- predict(fit, y.dim = 3, newdata = data.frame(dim1 = x_val), type = type)
    pred_23 <- predict(fit, y.dim = 3, newdata = data.frame(dim2 = x_val), type = type)

    # The preds shuld be dif when using different dims
    expect_false(identical(pred_12[1], pred_13[1]))
    expect_false(identical(pred_12[1], pred_23[1]))
    expect_false(identical(pred_13[1], pred_23[1]))
  }
})

# -------------------------------------------------------------------------------------------------------
