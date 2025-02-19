test_that("morphometric simulation correctly implements Ben's specifications", {
  # Bens parameters
  true_mus <- c(290, 125, 75)      # Mean dimensions
  true_sigmas <- c(45, 25, 15)     # Sd
  true_rhos <- c(0.75, 0.80, 0.85) # Corelations
  true_psis <- c(2.0, 1.5, 1.0)    # Scale parameters
  true_phis <- c(0.4, 0.5, 0.6)    # Error parameters

  # Parameter that create known relationship:
  # - Dims 2 and 3 should be isometric
  # - Dim 1 should have non-isometric relationships with others

  # Test setup
  set.seed(123)
  n.animals <- 50
  n.photos <- 5    # Fixed (per Ben's spec)
  m <- 3           # Num dims

  # Generate test data
  true_pars <- c(true_mus, true_sigmas, true_rhos, true_psis, true_phis)
  data <- sim.measurements(n.animals, rep(n.photos, n.animals), m, true_pars)

  # Fit model
  fit <- fit.morph(data)

  # 1. Test parameter recovery
  estimates <- summary(fit)
  expect_true(all(abs((estimates[1:3, "Estimate"] - true_mus)/true_mus) < 0.1),
              "Mean parameters should be recovered within 10%")

  # 2. Tst isometric relationships
  betas_23 <- summary(fit, type = "betas-pca", y.dim = 2, x.dim = 3)
  ratio_23 <- true_mus[2]/true_mus[3]
  expect_true(abs(betas_23[2,1] - ratio_23)/ratio_23 < 0.2,
              "Dimensions 2 and 3 should show isometric relationship")

  # 3. Test non-isometric relationships
  betas_12 <- summary(fit, type = "betas-pca", y.dim = 1, x.dim = 2)
  ratio_12 <- true_mus[1]/true_mus[2]
  expect_false(abs(betas_12[2,1] - ratio_12)/ratio_12 < 0.2,
               "Dimensions 1 and 2 should show non-isometric relationship")
})
