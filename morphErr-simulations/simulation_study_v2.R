# Load package
library(morphErr)

# Parameters frm manuscript
TRUE_PARAMS <- list(
  mus = c(315, 150, 100),      # Mean dimensions
  sigmas = c(25, 15, 10),     # Sd
  rhos = c(0.85, 0.80, 0.75), # Corelations between dimensions
  psis = c(10, 6, 4),    # Scale parameters
  phis = c(0.5, 0.4, 0.3)    # Error parameters
)

# Simulation settings
# Numbsimulations per sample size:
N_SIMS <- 1000
# Num picc per animal
N_PHOTOS <- 3
# Dif nums of animals to test
SAMPLE_SIZES <- c(15, 25, 50, 100)

cat("Starting simulation study...\n")

run_simulation_study <- function() {
  cat("Initialising simulation study\n")
  results <- list()

  # Loop through sample sizes
  for(n_animals in SAMPLE_SIZES) {
    cat(sprintf("\nRunning simulations for %d animals...\n", n_animals))

    # Run sims using sim.morph w/parallel processing
    sim_results <- sim.morph(
      n.sims = N_SIMS,
      n.animals = n_animals,
      n.photos = N_PHOTOS,
      mus = TRUE_PARAMS$mus,
      sigmas = TRUE_PARAMS$sigmas,
      rhos = TRUE_PARAMS$rhos,
      psis = TRUE_PARAMS$psis,
      phis = TRUE_PARAMS$phis,
      n.cores = 4,
      progressbar = TRUE
    )

    # Use eficient extraction for param estimates + ses
    parameter_estimates <- extract.sim.morph(sim_results, function(x) summary(x)[,1])
    parameter_ses <- extract.sim.morph(sim_results, function(x) summary(x)[,2])

    # Extract isometry tests efficiently
    isometry_tests <- extract.sim.morph(sim_results, function(x) {
      tests <- summary(x, type = "isometric-pca")
      c(
        # should be non-sig (isometric):
        tests["dim2 vs dim3", "P-value"],
        # should be sig (non-isometric):
        tests["dim1 vs dim2", "P-value"],
        # should be sig (non-isometric):
        tests["dim1 vs dim3", "P-value"]
      )
    })

    # Getdim relationships (LM + RMA)
    dim_relationships <- extract.sim.morph(sim_results, function(x) {
      # Get both interpretations for each pair
      c(
        # dims 2 vs 3 (should be isometric)
        summary(x, type = "betas", y.dim = 2, x.dim = 3)[2,1],
        summary(x, type = "betas-pca", y.dim = 2, x.dim = 3)[2,1],
        # dims 1 vs 2 (non-isometric)
        summary(x, type = "betas", y.dim = 1, x.dim = 2)[2,1],
        summary(x, type = "betas-pca", y.dim = 1, x.dim = 2)[2,1],
        # dims 1 vs 3 (non-isometric)
        summary(x, type = "betas", y.dim = 1, x.dim = 3)[2,1],
        summary(x, type = "betas-pca", y.dim = 1, x.dim = 3)[2,1]
      )
    })

    # Process true params
    true_pars <- c(TRUE_PARAMS$mus, TRUE_PARAMS$sigmas,
                   TRUE_PARAMS$rhos, TRUE_PARAMS$psis, TRUE_PARAMS$phis)

    # Calculate ci + coverage
    lower_ci <- parameter_estimates - 1.96 * parameter_ses
    upper_ci <- parameter_estimates + 1.96 * parameter_ses
    ci_coverage <- colMeans(true_pars >= lower_ci & true_pars <= upper_ci)

    # Store results for this sample size
    results[[as.character(n_animals)]] <- list(
      sample_size = n_animals,
      parameter_bias = colMeans(parameter_estimates) - true_pars,
      parameter_rmse = sqrt(colMeans((parameter_estimates - true_pars)^2)),
      ci_coverage = ci_coverage,
      isometry_power = c(
        mean(isometry_tests[,1] < 0.05),
        mean(isometry_tests[,2] < 0.05),
        mean(isometry_tests[,3] < 0.05)
      ),
      relationships = colMeans(dim_relationships)
    )

    # Print progress summary
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

    cat("\nMean Dimension Relationships (LM then RMA for each pair):\n")
    print(round(results[[as.character(n_animals)]]$relationships, 3))
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

# Run study if file is run directly
if (sys.nframe() == 0 || interactive()) {
  cat("Starting to run simulation study...\n")
  results <- run_simulation_study()
  cat("Simulation study complete!\n")
}
