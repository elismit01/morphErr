# Parameters frm manuscript
TRUE_PARAMS <- list(
  mus = c(290, 125, 75),      # Mean dimensions
  sigmas = c(45, 25, 15),     # Sd
  rhos = c(0.75, 0.80, 0.85), # Corelations between dimensions
  psis = c(2.0, 1.5, 1.0),    # Scale parameters
  phis = c(0.4, 0.5, 0.6)    # Error parameters
)

# Simulation settings
# Numbsimulations per sample size:
N_SIMS <- 1000
# Num picc per animal
N_PHOTOS <- 5
# Dif nums of animals to test
SAMPLE_SIZES <- c(20, 50, 100, 200)

cat("Starting simulation study...\n")

run_simulation_study <- function() {
  cat("Initialising simulation study\n")

  # Storage for results
  results <- list()

  # Loop through sample sizes
  for(n_animals in SAMPLE_SIZES) {
    cat(sprintf("\nRunning simulations for %d animals...\n", n_animals))

    # Storag for this sample size
    size_results <- list(
      parameter_estimates = matrix(NA, nrow = N_SIMS, ncol = 15),
      # Add storage for ses:
      parameter_ses = matrix(NA, nrow = N_SIMS, ncol = 15),
      # Track ci coverage:
      ci_coverage = matrix(FALSE, nrow = N_SIMS, ncol = 15),
      isometry_tests = matrix(NA, nrow = N_SIMS, ncol = 3)
    )

    # Run sims for this sample size
    for(i in 1:N_SIMS) {
      if(i %% 100 == 0) cat(sprintf("  Simulation %d of %d\n", i, N_SIMS))

      # Generate data
      true_pars <- c(TRUE_PARAMS$mus, TRUE_PARAMS$sigmas,
                     TRUE_PARAMS$rhos, TRUE_PARAMS$psis, TRUE_PARAMS$phis)
      data <- sim.measurements(n_animals, rep(N_PHOTOS, n_animals), m = 3, true_pars)

      # Fit model
      fit <- suppressWarnings(fit.morph(data))

      # Get sumary with estimates + ses
      fit_summary <- summary(fit)

      # Store parameter estimates + ses
      # Estimates:
      size_results$parameter_estimates[i,] <- fit_summary[,1]
      # Ses:
      size_results$parameter_ses[i,] <- fit_summary[,2]

      # Check if true vals fall within 95% cis
      lower_ci <- fit_summary[,1] - 1.96 * fit_summary[,2]
      upper_ci <- fit_summary[,1] + 1.96 * fit_summary[,2]
      size_results$ci_coverage[i,] <- true_pars >= lower_ci & true_pars <= upper_ci

      # Store isometry test results
      isometry_tests <- summary(fit, type = "isometric-pca")
      size_results$isometry_tests[i,] <- c(
        # should be non-sig (dims 2&3 are isometric):
        isometry_tests["dim2 vs dim3", "P-value"],
        # should be sig (non-isometric):
        isometry_tests["dim1 vs dim2", "P-value"],
        # should be sig (non-isometric):
        isometry_tests["dim1 vs dim3", "P-value"]
      )
    }

    # Calculate summary stats
    results[[as.character(n_animals)]] <- list(
      sample_size = n_animals,
      # Original summary stats
      parameter_bias = colMeans(size_results$parameter_estimates) - true_pars,
      parameter_rmse = sqrt(colMeans((size_results$parameter_estimates - true_pars)^2)),

      # New CI coverage summaries (prop of times true value in CI)
      ci_coverage = colMeans(size_results$ci_coverage),

      # Isometry test summaries
      isometry_power = c(
        # Type I error rate for 2vs3 (should be ~0.05):
        mean(size_results$isometry_tests[,1] < 0.05),
        # Power for 1vs2:
        mean(size_results$isometry_tests[,2] < 0.05),
        # Power for 1v3:
        mean(size_results$isometry_tests[,3] < 0.05)
      )
    )

    # Print summary for sample size
    cat(sprintf("\nResults for n = %d:\n", n_animals))
    cat("\nParameter Coverage Probabilities (should be close to 0.95):\n")
    names(results[[as.character(n_animals)]]$ci_coverage) <-
      c("mu1", "mu2", "mu3", "sigma1", "sigma2", "sigma3",
        "rho12", "rho13", "rho23", "psi1", "psi2", "psi3",
        "phi12", "phi13", "phi23")
    print(round(results[[as.character(n_animals)]]$ci_coverage, 3))

    cat("\nIsometry Test Results:\n")
    cat(sprintf("Type I Error Rate (2vs3): %.3f (should be ~0.05)\n",
                results[[as.character(n_animals)]]$isometry_power[1]))
    cat(sprintf("Power (1vs2): %.3f\n",
                results[[as.character(n_animals)]]$isometry_power[2]))
    cat(sprintf("Power (1vs3): %.3f\n",
                results[[as.character(n_animals)]]$isometry_power[3]))
  }

  # Add sim settings to results
  results$settings <- list(
    n_sims = N_SIMS,
    sample_sizes = SAMPLE_SIZES,
    n_photos = N_PHOTOS,
    true_params = TRUE_PARAMS
  )

  # Save
  save(results, file = "sim-study.RData")
  cat("\nSimulation complete! Results saved to sim-study.RData\n")

  return(results)
}

# Run the study if this file is run directly
if (sys.nframe() == 0 || interactive()) {
  cat("Starting to run simulation study...\n")
  results <- run_simulation_study()
  cat("Simulation study complete!\n")
}
