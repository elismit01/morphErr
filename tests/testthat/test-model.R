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
