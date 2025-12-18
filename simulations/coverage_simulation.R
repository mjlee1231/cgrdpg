# Full Coverage Probability Simulation
# Empirically verify that coverage rates match the nominal level (e.g., 95%)

library(cgrdpg)

# ------------------------------------------------------------------------------
# Data Generation Functions (same as quick test)
# ------------------------------------------------------------------------------
generate_latent_curve <- function(n) {
  t_vals <- (1:n) / n
  X <- matrix(0, n, 3)
  X[, 1] <- 0.15 * sin(2 * pi * t_vals) + 0.6
  X[, 2] <- 0.15 * cos(2 * pi * t_vals) + 0.6
  X[, 3] <- 0.15 * cos(4 * pi * t_vals)
  return(X)
}

generate_grdpg_data <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1) {
  X0 <- generate_latent_curve(n)
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))
  Y0 <- X0 %*% S
  P <- X0 %*% t(Y0)

  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0

  Z0 <- matrix(rnorm(p_cov * d, sd = 0.5), p_cov, d)
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 0.3), p_cov, n)

  list(A = A, B = B, X0 = X0, Y0 = Y0, Z0 = Z0, P = P, S = S)
}

compute_G_in <- function(i, X0, Y0, Z0, tau = 0.05) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)

  G_net <- matrix(0, d, d)
  for (j in 1:n) {
    G_net <- G_net + w[j] * outer(Y0[j, ], Y0[j, ])
  }

  G_cov <- crossprod(Z0)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

# ------------------------------------------------------------------------------
# Procrustes Alignment
# ------------------------------------------------------------------------------
procrustes_rotation <- function(X_hat, X0) {
  # Find rotation R that minimizes ||X0 - X_hat R||_F
  M <- t(X_hat) %*% X0
  svd_res <- svd(M)
  R <- svd_res$u %*% t(svd_res$v)
  return(R)
}

# ------------------------------------------------------------------------------
# Single Replicate Coverage Test
# ------------------------------------------------------------------------------
test_coverage_replicate <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1,
                                     vertices_to_test = NULL, alpha = 0.05,
                                     maxit = 5, verbose = FALSE) {
  # Generate data
  data <- generate_grdpg_data(n, p_cov, d, p_sig, q_sig)
  X0 <- data$X0
  Y0 <- data$Y0
  Z0 <- data$Z0

  # Fit model
  fit <- fit_grdpg_cov(data$A, data$B, d = d, p = p_sig, q = q_sig,
                       maxit = maxit, tol = 1e-3)

  if (!fit$converged && verbose) {
    warning(sprintf("Model did not converge in %d iterations", maxit))
  }

  # Align via Procrustes
  R <- procrustes_rotation(fit$X, X0)
  X_aligned <- fit$X %*% R

  # Test vertices
  if (is.null(vertices_to_test)) {
    vertices_to_test <- 1:n
  }

  n_vertices <- length(vertices_to_test)
  coverage_matrix <- matrix(FALSE, n_vertices, d)
  errors_matrix <- matrix(0, n_vertices, d)
  se_matrix <- matrix(0, n_vertices, d)

  z_crit <- qnorm(1 - alpha / 2)

  for (idx in 1:n_vertices) {
    i <- vertices_to_test[idx]

    # Compute covariance matrix
    G_in <- compute_G_in(i, X0, Y0, Z0, tau = fit$tau)

    # Check invertibility
    eigvals <- eigen(G_in, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigvals) < 1e-10) {
      if (verbose) {
        warning(sprintf("G_in is nearly singular for vertex %d", i))
      }
      next
    }

    Sigma_i <- solve(G_in)
    std_errors <- sqrt(diag(Sigma_i))
    se_matrix[idx, ] <- std_errors

    # Compute errors
    errors_matrix[idx, ] <- X_aligned[i, ] - X0[i, ]

    # Check coverage
    for (k in 1:d) {
      ci_lower <- X_aligned[i, k] - z_crit * std_errors[k]
      ci_upper <- X_aligned[i, k] + z_crit * std_errors[k]
      coverage_matrix[idx, k] <- (X0[i, k] >= ci_lower) && (X0[i, k] <= ci_upper)
    }
  }

  list(
    coverage = coverage_matrix,
    errors = errors_matrix,
    std_errors = se_matrix,
    converged = fit$converged,
    vertices_tested = vertices_to_test
  )
}

# ------------------------------------------------------------------------------
# Monte Carlo Simulation
# ------------------------------------------------------------------------------
run_coverage_simulation <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1,
                                     n_reps = 100, vertices_to_test = NULL,
                                     alpha = 0.05, seed = 123, verbose = TRUE) {
  set.seed(seed)

  # Start timing
  start_time <- Sys.time()

  if (verbose) {
    cat("=== Coverage Probability Simulation ===\n")
    cat(sprintf("Parameters: n=%d, p_cov=%d, d=%d, reps=%d, alpha=%.3f\n",
                n, p_cov, d, n_reps, alpha))
    cat(sprintf("Expected coverage: %.1f%%\n", 100 * (1 - alpha)))
    cat(sprintf("Start time: %s\n\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  }

  # If no vertices specified, test a sample
  if (is.null(vertices_to_test)) {
    # Test 10 evenly spaced vertices
    vertices_to_test <- round(seq(1, n, length.out = min(10, n)))
  }

  n_vertices <- length(vertices_to_test)

  # Storage for results
  coverage_count <- matrix(0, n_vertices, d)
  error_sum <- matrix(0, n_vertices, d)
  error_sq_sum <- matrix(0, n_vertices, d)
  se_sum <- matrix(0, n_vertices, d)
  n_converged <- 0
  replicate_times <- numeric(n_reps)

  # Run replicates
  for (rep in 1:n_reps) {
    rep_start <- Sys.time()

    if (verbose && rep %% 10 == 0) {
      elapsed <- as.numeric(difftime(rep_start, start_time, units = "secs"))
      eta <- (elapsed / rep) * (n_reps - rep)
      cat(sprintf("Replicate %d/%d (Elapsed: %.1fs, ETA: %.1fs)\n",
                  rep, n_reps, elapsed, eta))
    }

    result <- test_coverage_replicate(
      n, p_cov, d, p_sig, q_sig,
      vertices_to_test, alpha,
      verbose = FALSE
    )

    if (result$converged) {
      n_converged <- n_converged + 1
    }

    coverage_count <- coverage_count + result$coverage
    error_sum <- error_sum + result$errors
    error_sq_sum <- error_sq_sum + result$errors^2
    se_sum <- se_sum + result$std_errors

    # Record time for this replicate
    replicate_times[rep] <- as.numeric(difftime(Sys.time(), rep_start, units = "secs"))
  }

  # End timing
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Compute empirical coverage rates
  coverage_rates <- coverage_count / n_reps

  # Compute error statistics
  mean_errors <- error_sum / n_reps
  var_errors <- (error_sq_sum / n_reps) - mean_errors^2
  rmse <- sqrt(error_sq_sum / n_reps)
  mean_se <- se_sum / n_reps

  # Print results
  if (verbose) {
    cat("\n=== Results ===\n")
    cat(sprintf("Convergence rate: %.1f%%\n", 100 * n_converged / n_reps))
    cat("\nCoverage rates by coordinate:\n")
    for (k in 1:d) {
      cat(sprintf("  Coordinate %d: %.1f%% (SE: %.2f%%)\n",
                  k, 100 * mean(coverage_rates[, k]),
                  100 * sd(coverage_rates[, k]) / sqrt(n_vertices)))
    }
    cat(sprintf("\nOverall coverage: %.1f%%\n", 100 * mean(coverage_rates)))
    cat(sprintf("Expected: %.1f%%\n", 100 * (1 - alpha)))

    cat("\nAverage RMSE by coordinate:\n")
    for (k in 1:d) {
      cat(sprintf("  Coordinate %d: %.4f (avg SE: %.4f)\n",
                  k, mean(rmse[, k]), mean(mean_se[, k])))
    }

    cat("\n=== Timing ===\n")
    cat(sprintf("Total time: %.2f seconds (%.2f minutes)\n",
                total_time, total_time / 60))
    cat(sprintf("Average time per replicate: %.3f seconds\n",
                mean(replicate_times)))
    cat(sprintf("Median time per replicate: %.3f seconds\n",
                median(replicate_times)))
    cat(sprintf("Min/Max time per replicate: %.3f / %.3f seconds\n",
                min(replicate_times), max(replicate_times)))
    cat(sprintf("End time: %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
  }

  list(
    coverage_rates = coverage_rates,
    mean_errors = mean_errors,
    rmse = rmse,
    mean_se = mean_se,
    n_converged = n_converged,
    n_reps = n_reps,
    vertices_tested = vertices_to_test,
    params = list(n = n, p_cov = p_cov, d = d, p_sig = p_sig,
                  q_sig = q_sig, alpha = alpha),
    timing = list(
      total_time_seconds = total_time,
      mean_replicate_time = mean(replicate_times),
      median_replicate_time = median(replicate_times),
      min_replicate_time = min(replicate_times),
      max_replicate_time = max(replicate_times),
      all_replicate_times = replicate_times,
      start_time = start_time,
      end_time = end_time
    )
  )
}

# ------------------------------------------------------------------------------
# RUN SIMULATION
# ------------------------------------------------------------------------------
cat("Starting coverage simulation...\n\n")

# Main simulation
results <- run_coverage_simulation(
  n = 60,
  p_cov = 5,
  d = 3,
  p_sig = 2,
  q_sig = 1,
  n_reps = 200,
  vertices_to_test = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
  alpha = 0.05,
  seed = 123,
  verbose = TRUE
)

# Create output directory if it doesn't exist
if (!dir.exists("simulations")) {
  dir.create("simulations", recursive = TRUE)
  cat("\nCreated simulations directory\n")
}

# Save results
saveRDS(results, "simulations/coverage_results.rds")
cat("\nResults saved to simulations/coverage_results.rds\n")

# Optional: Test different sample sizes
cat("\n\n=== Testing different sample sizes ===\n")
sample_sizes <- c(40, 60, 80, 100)
size_results <- list()

for (n in sample_sizes) {
  cat(sprintf("\n--- n = %d ---\n", n))
  size_results[[as.character(n)]] <- run_coverage_simulation(
    n = n,
    p_cov = 5,
    d = 3,
    n_reps = 100,
    vertices_to_test = round(seq(1, n, length.out = 10)),
    alpha = 0.05,
    verbose = TRUE
  )
}

saveRDS(size_results, "simulations/coverage_by_size.rds")
cat("\nSize comparison saved to simulations/coverage_by_size.rds\n")
