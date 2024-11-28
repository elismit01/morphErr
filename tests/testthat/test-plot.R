# Create data once for all tests
create_test_data <- function() {
  set.seed(123)
  sim.measurements(
    n.animals = 10,
    n.photos = 3,
    m = 2,
    pars = c(100, 50,  # means
             10, 5,  # Sds
             0.7,  # correlation
             1, 0.5,  # measurement SDs
             0.3)  # measurement correlation
  )
}

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
  fit <- suppressWarnings(fit.morph(test_data))

    # Test diffferent plot types
  expect_silent(plot(fit, type = "data"))
  expect_silent(plot(fit, type = "ratio"))
})

test_that("plot.lme.morph handles line overlays", {
    # Create and fit test data
    test_data <- create_test_data()
    fit <- suppressWarnings(fit.morph(test_data))

    # Test different line types
    expect_warning(plot(fit, type = "data", line.type = "lm"))
    expect_warning(plot(fit, type = "data", line.type = "pca"))

    # Test with confidence intervals
    expect_warning(plot(fit, type = "data", line.type = "lm",
                        confints = TRUE))
})

test_that("plot.lme.morph handles adding to existing plots", {
  # Create + fit test data
  test_data <- create_test_data()
  fit <- suppressWarnings(fit.morph(test_data))

  # Test first plot
  expect_warning(
    plot(fit, type = "data", line.type = "lm"),
    "Could not add line to plot: Could not compute valid coefficients"
  )

  # Test ading to existing plot
  expect_warning(
    plot(fit, type = "data", line.type = "pca", add = TRUE, lty = 2),
    "Could not add line to plot: Could not compute valid coefficients"
  )
})

