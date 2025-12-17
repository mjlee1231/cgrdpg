# Coverage Simulation for GRDPG with Covariates
# Comprehensive simulation testing coverage rates with n=100 and n=500
# Parameters:
#   - vertices_to_test = n (test all vertices)
#   - p_cov = n/2 (number of covariates)

library(cgrdpg)

# Function to compute asymptotic variance-covariance matrix
compute_asymptotic_variance <- function(A, X, Z, B, sign_diag, tau = 0.05) {
  n <- nrow(X)
  d <- ncol(X)

  # Simplified asymptotic variance based on Fisher information
  Y <- X %*% sign_diag
  S <- X %*% t(Y)

  # Compute weights
  w <- dpsi(S, tau = tau)

  # Approximate variance for each vertex and dimension
  var_matrix <- matrix(0, n, d)
  for (i in 1:n) {
    # Network contribution to variance
    info_net <- sum(w[i, ] * rowSums(Y^2))

    # Covariate contribution
    info_cov <- sum(Z^2)

    # Total information
    info_total <- info_net + info_cov

    # Asymptotic variance (inverse information)
    if (info_total > 0) {
      var_matrix[i, ] <- 1 / info_total
    } else {
      var_matrix[i, ] <- Inf
    }
  }

  var_matrix
}

# Function to compute confidence intervals for all vertices
compute_all_vertex_ci <- function(X_hat, X_true, var_matrix, alpha = 0.05) {
  n <- nrow(X_hat)
  d <- ncol(X_hat)
  vertices_to_test <- n  # Test all vertices as specified

  z_crit <- qnorm(1 - alpha/2)

  coverage <- matrix(0, vertices_to_test, d)
  ci_lower <- matrix(0, vertices_to_test, d)
  ci_upper <- matrix(0, vertices_to_test, d)

  for (i in 1:vertices_to_test) {
    for (j in 1:d) {
      se <- sqrt(var_matrix[i, j])
      if (is.finite(se)) {
        ci_lower[i, j] <- X_hat[i, j] - z_crit * se
        ci_upper[i, j] <- X_hat[i, j] + z_crit * se
        coverage[i, j] <- (X_true[i, j] >= ci_lower[i, j]) &
                          (X_true[i, j] <= ci_upper[i, j])
      } else {
        coverage[i, j] <- NA
      }
    }
  }

  list(
    coverage = coverage,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    coverage_rate = mean(coverage, na.rm = TRUE),
    coverage_by_dim = colMeans(coverage, na.rm = TRUE),
    coverage_by_vertex = rowMeans(coverage, na.rm = TRUE)
  )
}

# Function to generate data and run one simulation
run_coverage_simulation <- function(n, d = 3, p_sig = 2, q_sig = 1,
                                   tau = 0.05, maxit = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Set p_cov = n/2 as specified
  p_cov <- n %/% 2
  vertices_to_test <- n  # Test all vertices

  # Generate true latent positions
  X_true <- matrix(rnorm(n * d, sd = 0.7), n, d)
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))

  # Generate adjacency matrix
  P <- X_true %*% S %*% t(X_true)
  P <- pmin(pmax(P, 1e-3), 1 - 1e-3)
  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0

  # Generate covariate matrix with Z
  Z_true <- matrix(rnorm(p_cov * d), p_cov, d)
  B <- Z_true %*% t(X_true) + matrix(rnorm(p_cov * n), p_cov, n)

  # Fit the model
  fit <- tryCatch({
    fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig,
                  tau = tau, maxit = maxit, tol = 1e-3)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(fit)) {
    return(list(success = FALSE))
  }

  # Align estimates to true values using Procrustes
  # X_hat W = X_true, solve for W
  XtX_inv <- tryCatch(solve(t(fit$X) %*% fit$X), error = function(e) NULL)
  if (is.null(XtX_inv)) {
    return(list(success = FALSE))
  }

  W <- XtX_inv %*% t(fit$X) %*% X_true
  X_aligned <- fit$X %*% W

  # Compute asymptotic variance
  var_matrix <- compute_asymptotic_variance(A, fit$X, fit$Z, B,
                                            diag(c(rep(1, p_sig), rep(-1, q_sig))),
                                            tau = tau)

  # Compute confidence intervals and coverage
  ci_result <- compute_all_vertex_ci(X_aligned, X_true, var_matrix, alpha = 0.05)

  list(
    success = TRUE,
    n = n,
    p_cov = p_cov,
    vertices_to_test = vertices_to_test,
    coverage_rate = ci_result$coverage_rate,
    coverage_by_dim = ci_result$coverage_by_dim,
    coverage_by_vertex = ci_result$coverage_by_vertex,
    converged = fit$converged,
    iters = fit$iters
  )
}

# Main function to run coverage simulation study
run_coverage_study <- function(n_vals = c(100, 500),
                               n_sims = 200,
                               parallel = FALSE) {
  results <- list()

  for (n in n_vals) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Running coverage simulation for n = %d\n", n))
    cat(sprintf("Parameters:\n"))
    cat(sprintf("  - vertices_to_test = %d (all vertices)\n", n))
    cat(sprintf("  - p_cov = %d (n/2)\n", n %/% 2))
    cat(sprintf("  - n_sims = %d\n", n_sims))
    cat(sprintf("========================================\n\n"))

    sim_results <- list()
    coverage_rates <- numeric(n_sims)
    coverage_by_dim <- matrix(0, n_sims, 3)
    converged_count <- 0
    success_count <- 0

    for (sim in 1:n_sims) {
      if (sim %% 20 == 0) {
        cat(sprintf("  Progress: %d/%d simulations completed\n", sim, n_sims))
      }

      result <- run_coverage_simulation(n = n, seed = sim * 1234 + n)

      if (result$success) {
        success_count <- success_count + 1
        coverage_rates[sim] <- result$coverage_rate
        coverage_by_dim[sim, ] <- result$coverage_by_dim
        if (result$converged) converged_count <- converged_count + 1
      } else {
        coverage_rates[sim] <- NA
        coverage_by_dim[sim, ] <- NA
      }

      sim_results[[sim]] <- result
    }

    # Compute summary statistics
    valid_coverage <- coverage_rates[!is.na(coverage_rates)]

    results[[paste0("n", n)]] <- list(
      n = n,
      p_cov = n %/% 2,
      vertices_to_test = n,
      n_sims = n_sims,
      success_rate = success_count / n_sims,
      convergence_rate = converged_count / success_count,
      coverage_rates = valid_coverage,
      mean_coverage = mean(valid_coverage),
      sd_coverage = sd(valid_coverage),
      median_coverage = median(valid_coverage),
      coverage_by_dim = coverage_by_dim,
      mean_coverage_by_dim = colMeans(coverage_by_dim, na.rm = TRUE),
      target_coverage = 0.95,
      sim_results = sim_results
    )

    # Print summary
    cat(sprintf("\nResults for n = %d:\n", n))
    cat(sprintf("  Success rate: %.2f%% (%d/%d)\n",
                100 * success_count / n_sims, success_count, n_sims))
    cat(sprintf("  Convergence rate: %.2f%% (%d/%d)\n",
                100 * converged_count / success_count, converged_count, success_count))
    cat(sprintf("  Mean coverage rate: %.4f (target: 0.95)\n", mean(valid_coverage)))
    cat(sprintf("  SD coverage rate: %.4f\n", sd(valid_coverage)))
    cat(sprintf("  Median coverage rate: %.4f\n", median(valid_coverage)))
    cat(sprintf("  Coverage by dimension:\n"))
    cat(sprintf("    Dim 1: %.4f\n", mean(coverage_by_dim[, 1], na.rm = TRUE)))
    cat(sprintf("    Dim 2: %.4f\n", mean(coverage_by_dim[, 2], na.rm = TRUE)))
    cat(sprintf("    Dim 3: %.4f\n\n", mean(coverage_by_dim[, 3], na.rm = TRUE)))
  }

  results
}

# Visualization function
plot_coverage_study <- function(results, output_file = "coverage_simulation.pdf") {
  pdf(output_file, width = 14, height = 10)

  par(mfrow = c(2, 3))

  n_vals <- sapply(results, function(x) x$n)

  # Plot 1: Overall coverage rate comparison
  coverage_means <- sapply(results, function(x) x$mean_coverage)
  coverage_sds <- sapply(results, function(x) x$sd_coverage)

  plot(n_vals, coverage_means,
       type = "b", pch = 19, col = "blue", cex = 1.5,
       ylim = c(0.8, 1.0),
       xlab = "Sample size (n)",
       ylab = "Coverage rate",
       main = "Overall Coverage Rate vs Sample Size",
       cex.lab = 1.2, cex.main = 1.3)
  arrows(n_vals, coverage_means - 2*coverage_sds,
         n_vals, coverage_means + 2*coverage_sds,
         angle = 90, code = 3, length = 0.1, col = "blue", lwd = 2)
  abline(h = 0.95, col = "red", lty = 2, lwd = 2)
  legend("bottomright",
         legend = c("Observed ± 2 SD", "Target (95%)"),
         col = c("blue", "red"), lty = c(1, 2), pch = c(19, NA), cex = 1.1)

  # Plot 2: Distribution comparison
  boxplot(lapply(results, function(x) x$coverage_rates),
          names = paste0("n=", n_vals),
          col = c("lightblue", "lightgreen"),
          ylab = "Coverage rate",
          main = "Distribution of Coverage Rates",
          cex.lab = 1.2, cex.main = 1.3)
  abline(h = 0.95, col = "red", lty = 2, lwd = 2)

  # Plot 3: Coverage by dimension for n=100
  if ("n100" %in% names(results)) {
    cov_dim <- results$n100$mean_coverage_by_dim
    barplot(cov_dim,
            names.arg = c("Dim 1", "Dim 2", "Dim 3"),
            col = "coral",
            ylim = c(0, 1),
            main = sprintf("Coverage by Dimension (n=100, p_cov=%d)",
                          results$n100$p_cov),
            ylab = "Coverage rate",
            cex.lab = 1.2, cex.main = 1.3)
    abline(h = 0.95, col = "red", lty = 2, lwd = 2)
  }

  # Plot 4: Coverage by dimension for n=500
  if ("n500" %in% names(results)) {
    cov_dim <- results$n500$mean_coverage_by_dim
    barplot(cov_dim,
            names.arg = c("Dim 1", "Dim 2", "Dim 3"),
            col = "lightgreen",
            ylim = c(0, 1),
            main = sprintf("Coverage by Dimension (n=500, p_cov=%d)",
                          results$n500$p_cov),
            ylab = "Coverage rate",
            cex.lab = 1.2, cex.main = 1.3)
    abline(h = 0.95, col = "red", lty = 2, lwd = 2)
  }

  # Plot 5: Success and convergence rates
  success_rates <- sapply(results, function(x) x$success_rate)
  conv_rates <- sapply(results, function(x) x$convergence_rate)

  barplot(rbind(success_rates, conv_rates),
          beside = TRUE,
          names.arg = paste0("n=", n_vals),
          col = c("skyblue", "orange"),
          ylim = c(0, 1),
          main = "Success and Convergence Rates",
          ylab = "Rate",
          cex.lab = 1.2, cex.main = 1.3)
  legend("bottomright",
         legend = c("Success", "Convergence"),
         fill = c("skyblue", "orange"),
         cex = 1.1)

  # Plot 6: Histogram of coverage for n=500
  if ("n500" %in% names(results)) {
    hist(results$n500$coverage_rates,
         breaks = 30,
         col = "lightgreen",
         main = sprintf("Coverage Rate Distribution (n=500)"),
         xlab = "Coverage rate",
         cex.lab = 1.2, cex.main = 1.3)
    abline(v = 0.95, col = "red", lty = 2, lwd = 2)
    abline(v = mean(results$n500$coverage_rates), col = "blue", lty = 2, lwd = 2)
    legend("topleft",
           legend = c("Target", "Observed mean"),
           col = c("red", "blue"),
           lty = 2,
           lwd = 2,
           cex = 1.1)
  }

  par(mfrow = c(1, 1))
  dev.off()

  cat(sprintf("\nPlots saved to: %s\n", output_file))
}

# Main execution
if (interactive() || !interactive()) {
  cat("================================================\n")
  cat("Coverage Simulation for GRDPG with Covariates\n")
  cat("================================================\n")
  cat("Testing with n = 100 and n = 500\n")
  cat("Parameters:\n")
  cat("  - vertices_to_test = n (all vertices)\n")
  cat("  - p_cov = n/2\n")
  cat("  - Target coverage: 95%\n")
  cat("================================================\n")

  # Run the coverage study
  results <- run_coverage_study(n_vals = c(100, 500), n_sims = 200)

  # Save results
  saveRDS(results, file = "coverage_simulation_results.rds")
  cat("\nResults saved to: coverage_simulation_results.rds\n")

  # Create visualizations
  plot_coverage_study(results, output_file = "coverage_simulation.pdf")

  # Print final summary
  cat("\n================================================\n")
  cat("FINAL SUMMARY\n")
  cat("================================================\n")
  for (i in seq_along(results)) {
    r <- results[[i]]
    cat(sprintf("\nn = %d (p_cov = %d, vertices_to_test = %d):\n",
                r$n, r$p_cov, r$vertices_to_test))
    cat(sprintf("  Mean coverage: %.4f ± %.4f (target: %.4f)\n",
                r$mean_coverage, r$sd_coverage, r$target_coverage))
    cat(sprintf("  Median coverage: %.4f\n", r$median_coverage))
    cat(sprintf("  Success rate: %.2f%%\n", 100 * r$success_rate))
    cat(sprintf("  Convergence rate: %.2f%%\n", 100 * r$convergence_rate))
  }
  cat("================================================\n")
}
