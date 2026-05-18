# Quick test of the simulation script with reduced repetitions
# This will run just 2 repetitions to verify everything works

library(cgrdpg)

# Source the main simulation functions
source("simulations_latent_curve.R", echo = FALSE)

cat("\n=== QUICK TEST RUN ===\n")
cat("Testing with n=100, p_cov=50, 2 repetitions\n\n")

# Run a quick test with just 2 repetitions
test_result <- run_simulation(n = 100, p_cov = 50, n_reps = 2, seed_start = 9999)

if (!is.na(test_result$mean_sse)) {
  cat("\n✓ Test completed successfully!\n")
  cat("The full simulation script is ready to run.\n\n")
} else {
  cat("\n✗ Test failed. Please check the error messages above.\n\n")
}
