# GRDPG Simulation Experiments with Latent Curve Parameterization
# =================================================================
#
# This script runs simulation experiments for a rank-three generic GRDPG
# with signature (2,1), where latent positions are drawn from a 3D curve.
#
# Latent curve: [0.15*sin(2*pi*t) + 0.6, 0.15*cos(2*pi*t) + 0.6, 0.15*cos(4*pi*t)]^T
# where t in [0,1], with n positions at t = i/n for i in [n]
#
# Scenarios: n = 100, 500 with p_cov = 0.5n
# Repetitions: 100 for each scenario

library(cgrdpg)

# Helper function: Generate latent positions from the parametric curve
generate_latent_curve <- function(n) {
  # t = i/n for i in [n], so t ranges from 1/n to 1
  t_vals <- (1:n) / n

  X <- matrix(0, n, 3)
  X[, 1] <- 0.15 * sin(2 * pi * t_vals) + 0.6
  X[, 2] <- 0.15 * cos(2 * pi * t_vals) + 0.6
  X[, 3] <- 0.15 * cos(4 * pi * t_vals)

  return(X)
}

# Helper function: Procrustes alignment for computing SSE
# Aligns estimated X_hat to true X0 via orthogonal transformation
procrustes_align <- function(X_hat, X0) {
  # Find Q such that ||X_hat Q - X0||_F is minimized
  # Solution: Q from SVD of t(X_hat) %*% X0
  svd_result <- svd(t(X_hat) %*% X0)
  Q <- svd_result$u %*% t(svd_result$v)

  X_aligned <- X_hat %*% Q
  return(X_aligned)
}

# Helper function: Compute SSE after Procrustes alignment
compute_sse <- function(X_hat, X0) {
  X_aligned <- procrustes_align(X_hat, X0)
  sse <- sum((X_aligned - X0)^2)
  return(sse)
}

# Main simulation function for a single scenario
run_simulation <- function(n, p_cov, n_reps = 100, d = 3, p_sig = 2, q_sig = 1,
                          maxit = 5, tol = 1e-2, tau = 0.05, seed_start = 1000) {

  cat(sprintf("\n=== Running simulation for n=%d, p_cov=%d ===\n", n, p_cov))
  cat(sprintf("Repetitions: %d\n", n_reps))
  cat(sprintf("Signature: (%d,%d), Embedding dimension: %d\n\n", p_sig, q_sig, d))

  # Storage for results
  sse_values <- numeric(n_reps)
  time_values <- numeric(n_reps)

  # Signature matrix
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))

  # Run repetitions
  for (rep in 1:n_reps) {
    if (rep %% 10 == 0) {
      cat(sprintf("  Repetition %d/%d...\n", rep, n_reps))
    }

    set.seed(seed_start + rep)

    # Generate true latent positions from the curve
    X0 <- generate_latent_curve(n)

    # Generate probability matrix
    P <- X0 %*% S %*% t(X0)
    # Clip probabilities to valid range [1e-3, 1-1e-3]
    P <- pmin(pmax(P, 1e-3), 1 - 1e-3)

    # Sample adjacency matrix (symmetric, no self-loops)
    A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
    A[lower.tri(A)] <- t(A)[lower.tri(A)]
    diag(A) <- 0

    # Generate covariate loadings Z0
    Z0 <- matrix(rnorm(p_cov * d, mean = 0, sd = 0.5), p_cov, d)

    # Generate covariate matrix B = Z0 X0^T + noise
    B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, mean = 0, sd = 1), p_cov, n)

    # Fit the model and time it
    start_time <- Sys.time()

    fit <- tryCatch({
      fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig,
                    maxit = maxit, tol = tol, tau = tau)
    }, error = function(e) {
      cat(sprintf("  ERROR in repetition %d: %s\n", rep, e$message))
      return(NULL)
    })

    end_time <- Sys.time()
    time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Compute SSE if fitting succeeded
    if (!is.null(fit) && !is.null(fit$X)) {
      sse <- compute_sse(fit$X, X0)
      sse_values[rep] <- sse
      time_values[rep] <- time_elapsed
    } else {
      # Mark as NA if fitting failed
      sse_values[rep] <- NA
      time_values[rep] <- NA
    }
  }

  # Compute summary statistics (removing NAs if any)
  valid_indices <- !is.na(sse_values)
  n_valid <- sum(valid_indices)

  if (n_valid == 0) {
    cat("\nWARNING: All repetitions failed!\n")
    return(list(
      n = n,
      p_cov = p_cov,
      n_reps = n_reps,
      n_valid = 0,
      mean_sse = NA,
      sd_sse = NA,
      mean_time = NA
    ))
  }

  mean_sse <- mean(sse_values[valid_indices])
  sd_sse <- sd(sse_values[valid_indices])
  mean_time <- mean(time_values[valid_indices])

  # Print results
  cat(sprintf("\n--- Results for n=%d, p_cov=%d ---\n", n, p_cov))
  cat(sprintf("Valid repetitions: %d/%d\n", n_valid, n_reps))
  cat(sprintf("Mean SSE:          %.6f\n", mean_sse))
  cat(sprintf("Std Dev SSE:       %.6f\n", sd_sse))
  cat(sprintf("Mean Runtime:      %.4f seconds\n", mean_time))
  cat(sprintf("Total Runtime:     %.2f seconds\n", sum(time_values[valid_indices])))
  cat("=====================================\n\n")

  return(list(
    n = n,
    p_cov = p_cov,
    n_reps = n_reps,
    n_valid = n_valid,
    mean_sse = mean_sse,
    sd_sse = sd_sse,
    mean_time = mean_time,
    sse_values = sse_values[valid_indices],
    time_values = time_values[valid_indices]
  ))
}

# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("===============================================================\n")
cat("  GRDPG Simulation with Latent Curve Parameterization\n")
cat("===============================================================\n")
cat("\nLatent curve: [0.15*sin(2*pi*t) + 0.6, \n")
cat("               0.15*cos(2*pi*t) + 0.6, \n")
cat("               0.15*cos(4*pi*t)]^T\n")
cat("\nSignature: (2, 1)\n")
cat("Repetitions per scenario: 100\n")
cat("===============================================================\n")

# Record overall start time
overall_start <- Sys.time()

# Scenario 1: n = 100, p_cov = 50
results_100 <- run_simulation(n = 100, p_cov = 50, n_reps = 100)

# Scenario 2: n = 500, p_cov = 250
results_500 <- run_simulation(n = 500, p_cov = 250, n_reps = 100)

# Overall timing
overall_end <- Sys.time()
total_time <- as.numeric(difftime(overall_end, overall_start, units = "mins"))

# =============================================================================
# Summary Report
# =============================================================================

cat("\n\n")
cat("===============================================================\n")
cat("                    FINAL SUMMARY REPORT\n")
cat("===============================================================\n\n")

summary_df <- data.frame(
  Scenario = c("n=100, p_cov=50", "n=500, p_cov=250"),
  N_Valid = c(results_100$n_valid, results_500$n_valid),
  Mean_SSE = c(results_100$mean_sse, results_500$mean_sse),
  SD_SSE = c(results_100$sd_sse, results_500$sd_sse),
  Mean_Runtime_sec = c(results_100$mean_time, results_500$mean_time)
)

print(summary_df)

cat("\n")
cat(sprintf("Total simulation time: %.2f minutes (%.2f hours)\n",
            total_time, total_time / 60))
cat("===============================================================\n\n")

# Save results to RDS file
results <- list(
  scenario_1 = results_100,
  scenario_2 = results_500,
  summary = summary_df,
  total_time_minutes = total_time,
  timestamp = Sys.time()
)

output_file <- "simulation_results_latent_curve.rds"
saveRDS(results, file = output_file)
cat(sprintf("Results saved to: %s\n\n", output_file))

# Also save as CSV for easy viewing
csv_file <- "simulation_summary_latent_curve.csv"
write.csv(summary_df, file = csv_file, row.names = FALSE)
cat(sprintf("Summary table saved to: %s\n\n", csv_file))
