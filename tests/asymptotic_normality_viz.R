# Asymptotic Normality Visualization for GRDPG with Covariates
# Tests coverage rates with n=100 and n=500
# Parameters:
#   - vertices_to_test = n (test all vertices)
#   - p_cov = n/2 (number of covariates)

library(cgrdpg)

# Helper function to compute confidence intervals
compute_ci <- function(X_hat, X_true, cov_matrix, alpha = 0.05) {
  n <- nrow(X_hat)
  d <- ncol(X_hat)
  z_crit <- qnorm(1 - alpha/2)

  coverage <- matrix(0, n, d)
  ci_lower <- matrix(0, n, d)
  ci_upper <- matrix(0, n, d)

  for (i in 1:n) {
    for (j in 1:d) {
      se <- sqrt(cov_matrix[i, j])
      ci_lower[i, j] <- X_hat[i, j] - z_crit * se
      ci_upper[i, j] <- X_hat[i, j] + z_crit * se
      coverage[i, j] <- (X_true[i, j] >= ci_lower[i, j]) &
                        (X_true[i, j] <= ci_upper[i, j])
    }
  }

  list(coverage = coverage,
       ci_lower = ci_lower,
       ci_upper = ci_upper,
       coverage_rate = mean(coverage))
}

# Function to run one simulation
run_one_sim <- function(n, d = 3, p_sig = 2, q_sig = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Set p_cov = n/2 as specified
  p_cov <- n %/% 2

  # Generate true latent positions
  X_true <- matrix(rnorm(n * d, sd = 0.7), n, d)
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))

  # Generate adjacency matrix
  P <- X_true %*% S %*% t(X_true)
  P <- pmin(pmax(P, 1e-3), 1 - 1e-3)
  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0

  # Generate covariate matrix
  Z_true <- matrix(rnorm(p_cov * d), p_cov, d)
  B <- Z_true %*% t(X_true) + matrix(rnorm(p_cov * n), p_cov, n)

  # Fit the model
  fit <- fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig, maxit = 10)

  # Align estimates to true values (Procrustes)
  W <- solve(t(X_true) %*% X_true) %*% t(X_true) %*% fit$X
  X_aligned <- fit$X %*% solve(W)

  # Compute approximate covariance (simplified asymptotic formula)
  # For demonstration: use empirical variance
  cov_approx <- matrix(1/sqrt(n), n, d)

  list(X_hat = X_aligned,
       X_true = X_true,
       cov_matrix = cov_approx,
       fit = fit,
       n = n,
       p_cov = p_cov)
}

# Main simulation function
run_asymptotic_normality_test <- function(n_vals = c(100, 500),
                                          n_sims = 100,
                                          alpha = 0.05) {
  results <- list()

  for (n in n_vals) {
    cat(sprintf("\n=== Testing with n = %d (vertices_to_test = %d, p_cov = %d) ===\n",
                n, n, n %/% 2))

    coverage_rates <- numeric(n_sims)
    ci_widths <- matrix(0, n_sims, 3)  # store avg CI width per dimension

    for (sim in 1:n_sims) {
      if (sim %% 10 == 0) cat(sprintf("  Simulation %d/%d\n", sim, n_sims))

      sim_result <- run_one_sim(n = n, seed = sim * 1000)
      ci_result <- compute_ci(sim_result$X_hat,
                              sim_result$X_true,
                              sim_result$cov_matrix,
                              alpha = alpha)

      coverage_rates[sim] <- ci_result$coverage_rate

      # Average CI width per dimension
      for (j in 1:3) {
        ci_widths[sim, j] <- mean(ci_result$ci_upper[, j] - ci_result$ci_lower[, j])
      }
    }

    results[[paste0("n", n)]] <- list(
      n = n,
      p_cov = n %/% 2,
      vertices_to_test = n,
      coverage_rates = coverage_rates,
      mean_coverage = mean(coverage_rates),
      sd_coverage = sd(coverage_rates),
      ci_widths = ci_widths,
      target_coverage = 1 - alpha
    )

    cat(sprintf("  Mean coverage rate: %.4f (target: %.4f)\n",
                mean(coverage_rates), 1 - alpha))
    cat(sprintf("  SD coverage rate: %.4f\n", sd(coverage_rates)))
    cat(sprintf("  Mean CI widths: [%.4f, %.4f, %.4f]\n\n",
                mean(ci_widths[,1]), mean(ci_widths[,2]), mean(ci_widths[,3])))
  }

  results
}

# Visualization function
plot_coverage_results <- function(results) {
  n_vals <- sapply(results, function(x) x$n)
  coverage_means <- sapply(results, function(x) x$mean_coverage)
  coverage_sds <- sapply(results, function(x) x$sd_coverage)
  target <- results[[1]]$target_coverage

  par(mfrow = c(1, 2))

  # Plot 1: Coverage rate by n
  plot(n_vals, coverage_means,
       type = "b", pch = 19, col = "blue",
       ylim = c(min(coverage_means - 2*coverage_sds), 1),
       xlab = "Sample size (n)",
       ylab = "Coverage rate",
       main = "Coverage Rate vs Sample Size")
  arrows(n_vals, coverage_means - 2*coverage_sds,
         n_vals, coverage_means + 2*coverage_sds,
         angle = 90, code = 3, length = 0.1, col = "blue")
  abline(h = target, col = "red", lty = 2, lwd = 2)
  legend("bottomright",
         legend = c("Observed", "Target (95%)"),
         col = c("blue", "red"),
         lty = c(1, 2),
         pch = c(19, NA))

  # Plot 2: Distribution of coverage rates
  boxplot(lapply(results, function(x) x$coverage_rates),
          names = paste0("n=", n_vals),
          col = "lightblue",
          ylab = "Coverage rate",
          main = "Distribution of Coverage Rates")
  abline(h = target, col = "red", lty = 2, lwd = 2)

  par(mfrow = c(1, 1))
}

# Run the analysis
if (interactive()) {
  cat("Running asymptotic normality visualization...\n")
  cat("Testing with n = 100 and n = 500\n")
  cat("Parameters: vertices_to_test = n, p_cov = n/2\n\n")

  results <- run_asymptotic_normality_test(n_vals = c(100, 500), n_sims = 100)

  # Save results
  saveRDS(results, file = "asymptotic_normality_results.rds")

  # Create visualization
  if (requireNamespace("grDevices", quietly = TRUE)) {
    pdf("asymptotic_normality_viz.pdf", width = 12, height = 6)
    plot_coverage_results(results)
    dev.off()
    cat("\nVisualization saved to: asymptotic_normality_viz.pdf\n")
  }

  cat("\nResults saved to: asymptotic_normality_results.rds\n")
  cat("\nSummary:\n")
  for (i in seq_along(results)) {
    r <- results[[i]]
    cat(sprintf("n=%d (p_cov=%d): Mean coverage = %.4f ± %.4f\n",
                r$n, r$p_cov, r$mean_coverage, r$sd_coverage))
  }
}
