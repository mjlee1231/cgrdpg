# Quick Coverage Test for GRDPG with Covariates
# Fast test script for coverage rates with n=100 and n=500
# Parameters:
#   - vertices_to_test = n (test all vertices)
#   - p_cov = n/2 (number of covariates)

library(cgrdpg)

# Quick helper to compute simple confidence intervals
quick_compute_ci <- function(X_hat, X_true, n, alpha = 0.05) {
  d <- ncol(X_hat)
  z_crit <- qnorm(1 - alpha/2)

  # Simple asymptotic SE approximation: 1/sqrt(n)
  se <- 1/sqrt(n)

  coverage <- matrix(0, nrow(X_hat), d)
  for (i in 1:nrow(X_hat)) {
    for (j in 1:d) {
      ci_lower <- X_hat[i, j] - z_crit * se
      ci_upper <- X_hat[i, j] + z_crit * se
      coverage[i, j] <- (X_true[i, j] >= ci_lower) & (X_true[i, j] <= ci_upper)
    }
  }

  list(coverage_rate = mean(coverage),
       coverage_by_dim = colMeans(coverage))
}

# Quick simulation function
quick_sim <- function(n, d = 3, p_sig = 2, q_sig = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Set parameters as specified
  p_cov <- n %/% 2           # p_cov = n/2
  vertices_to_test <- n      # Test all vertices

  # Generate data
  X_true <- matrix(rnorm(n * d, sd = 0.7), n, d)
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))

  P <- X_true %*% S %*% t(X_true)
  P <- pmin(pmax(P, 1e-3), 1 - 1e-3)
  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0

  Z_true <- matrix(rnorm(p_cov * d), p_cov, d)
  B <- Z_true %*% t(X_true) + matrix(rnorm(p_cov * n), p_cov, n)

  # Fit model
  fit <- tryCatch({
    fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig, maxit = 5, tol = 1e-3)
  }, error = function(e) NULL)

  if (is.null(fit)) return(list(success = FALSE))

  # Align using Procrustes
  W <- tryCatch({
    solve(t(fit$X) %*% fit$X) %*% t(fit$X) %*% X_true
  }, error = function(e) NULL)

  if (is.null(W)) return(list(success = FALSE))

  X_aligned <- fit$X %*% W

  # Quick coverage check
  ci_result <- quick_compute_ci(X_aligned, X_true, n)

  list(success = TRUE,
       n = n,
       p_cov = p_cov,
       vertices_to_test = vertices_to_test,
       coverage_rate = ci_result$coverage_rate,
       coverage_by_dim = ci_result$coverage_by_dim,
       converged = fit$converged)
}

# Quick test function
quick_coverage_test <- function(n_vals = c(100, 500), n_sims = 50) {
  cat("\n========================================\n")
  cat("QUICK COVERAGE TEST\n")
  cat("========================================\n")
  cat(sprintf("Testing with n = %s\n", paste(n_vals, collapse = ", ")))
  cat(sprintf("Parameters:\n"))
  cat(sprintf("  - vertices_to_test = n (all vertices)\n"))
  cat(sprintf("  - p_cov = n/2\n"))
  cat(sprintf("  - n_sims = %d (quick test)\n", n_sims))
  cat("========================================\n\n")

  results <- list()

  for (n in n_vals) {
    cat(sprintf("Testing n = %d (p_cov = %d)...\n", n, n %/% 2))

    coverage_rates <- numeric(n_sims)
    coverage_by_dim <- matrix(0, n_sims, 3)
    success_count <- 0
    converged_count <- 0

    for (sim in 1:n_sims) {
      if (sim %% 10 == 0) cat(sprintf("  %d/%d ", sim, n_sims))

      result <- quick_sim(n = n, seed = sim * 999)

      if (result$success) {
        success_count <- success_count + 1
        coverage_rates[sim] <- result$coverage_rate
        coverage_by_dim[sim, ] <- result$coverage_by_dim
        if (result$converged) converged_count <- converged_count + 1
      } else {
        coverage_rates[sim] <- NA
      }
    }
    cat("\n")

    valid_coverage <- coverage_rates[!is.na(coverage_rates)]

    results[[paste0("n", n)]] <- list(
      n = n,
      p_cov = n %/% 2,
      vertices_to_test = n,
      success_rate = success_count / n_sims,
      convergence_rate = if(success_count > 0) converged_count / success_count else 0,
      mean_coverage = mean(valid_coverage),
      sd_coverage = sd(valid_coverage),
      coverage_by_dim = colMeans(coverage_by_dim, na.rm = TRUE)
    )

    # Print quick summary
    cat(sprintf("\n  Results (n=%d, p_cov=%d, vertices_to_test=%d):\n",
                n, n %/% 2, n))
    cat(sprintf("    Success: %.1f%% | Converged: %.1f%%\n",
                100 * success_count / n_sims,
                100 * converged_count / success_count))
    cat(sprintf("    Coverage: %.4f ± %.4f (target: 0.95)\n",
                mean(valid_coverage), sd(valid_coverage)))
    cat(sprintf("    By dimension: [%.4f, %.4f, %.4f]\n\n",
                colMeans(coverage_by_dim, na.rm = TRUE)[1],
                colMeans(coverage_by_dim, na.rm = TRUE)[2],
                colMeans(coverage_by_dim, na.rm = TRUE)[3]))
  }

  results
}

# Simple visualization
plot_quick_results <- function(results) {
  n_vals <- sapply(results, function(x) x$n)
  coverage_means <- sapply(results, function(x) x$mean_coverage)
  coverage_sds <- sapply(results, function(x) x$sd_coverage)

  par(mfrow = c(1, 2))

  # Plot 1: Coverage by n
  plot(n_vals, coverage_means,
       type = "b", pch = 19, col = "darkgreen", cex = 2,
       ylim = c(0.85, 1),
       xlab = "Sample size (n)",
       ylab = "Coverage rate",
       main = "Quick Coverage Test Results",
       cex.lab = 1.2, cex.main = 1.3)
  arrows(n_vals, coverage_means - 2*coverage_sds,
         n_vals, coverage_means + 2*coverage_sds,
         angle = 90, code = 3, length = 0.1, col = "darkgreen", lwd = 2)
  abline(h = 0.95, col = "red", lty = 2, lwd = 2)
  text(n_vals, coverage_means + 0.02,
       labels = sprintf("%.3f", coverage_means),
       pos = 3, cex = 1.1)
  legend("bottomright",
         legend = c("Observed", "Target (95%)"),
         col = c("darkgreen", "red"),
         lty = c(1, 2),
         pch = c(19, NA))

  # Plot 2: Coverage by dimension
  dims_data <- t(sapply(results, function(x) x$coverage_by_dim))
  barplot(t(dims_data),
          beside = TRUE,
          names.arg = paste0("n=", n_vals),
          col = c("coral", "skyblue", "lightgreen"),
          ylim = c(0, 1),
          main = "Coverage by Dimension",
          ylab = "Coverage rate",
          cex.lab = 1.2, cex.main = 1.3)
  abline(h = 0.95, col = "red", lty = 2, lwd = 2)
  legend("bottomright",
         legend = c("Dim 1", "Dim 2", "Dim 3", "Target"),
         fill = c("coral", "skyblue", "lightgreen", NA),
         border = c("black", "black", "black", NA),
         col = c(NA, NA, NA, "red"),
         lty = c(NA, NA, NA, 2),
         lwd = c(NA, NA, NA, 2))

  par(mfrow = c(1, 1))
}

# Main execution
cat("\n")
cat("================================================\n")
cat("  QUICK COVERAGE TEST FOR GRDPG WITH COVARIATES\n")
cat("================================================\n")

# Run quick test
start_time <- Sys.time()
results <- quick_coverage_test(n_vals = c(100, 500), n_sims = 50)
end_time <- Sys.time()

# Save results
saveRDS(results, file = "test_coverage_quick_results.rds")
cat(sprintf("\nResults saved to: test_coverage_quick_results.rds\n"))

# Create visualization
if (requireNamespace("grDevices", quietly = TRUE)) {
  pdf("test_coverage_quick.pdf", width = 12, height = 6)
  plot_quick_results(results)
  dev.off()
  cat("Visualization saved to: test_coverage_quick.pdf\n")
}

# Print summary
cat("\n================================================\n")
cat("QUICK TEST SUMMARY\n")
cat("================================================\n")
cat(sprintf("Time elapsed: %.2f seconds\n\n", as.numeric(end_time - start_time, units = "secs")))

for (i in seq_along(results)) {
  r <- results[[i]]
  cat(sprintf("n = %d (p_cov = %d, vertices_to_test = %d):\n",
              r$n, r$p_cov, r$vertices_to_test))
  cat(sprintf("  Coverage: %.4f ± %.4f (target: 0.95)\n",
              r$mean_coverage, r$sd_coverage))
  cat(sprintf("  Success rate: %.1f%%\n", 100 * r$success_rate))
  cat(sprintf("  Convergence rate: %.1f%%\n", 100 * r$convergence_rate))
  cat(sprintf("  Dimension coverage: [%.4f, %.4f, %.4f]\n\n",
              r$coverage_by_dim[1], r$coverage_by_dim[2], r$coverage_by_dim[3]))
}
cat("================================================\n")
cat("\nQuick test complete!\n")
cat("For comprehensive results, run coverage_simulation.R\n")
cat("================================================\n\n")
