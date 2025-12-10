# Quick test for latent position coverage probability
# Before running full simulations, verify the setup is viable

library(cgrdpg)

# ------------------------------------------------------------------------------
# 1. Generate true latent positions using the curve
# ------------------------------------------------------------------------------
generate_latent_curve <- function(n) {
  # t = i/n for i in [n], so t ranges from 1/n to 1
  t_vals <- (1:n) / n

  X <- matrix(0, n, 3)
  X[, 1] <- 0.15 * sin(2 * pi * t_vals) + 0.6
  X[, 2] <- 0.15 * cos(2 * pi * t_vals) + 0.6
  X[, 3] <- 0.15 * cos(4 * pi * t_vals)

  return(X)
}

# ------------------------------------------------------------------------------
# 2. Generate network and covariates from GRDPG model
# ------------------------------------------------------------------------------
generate_grdpg_data <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1) {
  # True latent positions
  X0 <- generate_latent_curve(n)

  # Signature matrix
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))

  # Compute probability matrix
  Y0 <- X0 %*% S
  P <- X0 %*% t(Y0)

  # Clip probabilities to avoid numerical issues
  P <- pmin(pmax(P, 1e-3), 1 - 1e-3)

  # Generate adjacency matrix
  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]  # Symmetrize
  diag(A) <- 0  # No self-loops

  # Generate covariate loadings Z
  Z0 <- matrix(rnorm(p_cov * d, sd = 0.5), p_cov, d)

  # Generate covariates: B = Z0 X0^T + noise
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 0.3), p_cov, n)

  list(A = A, B = B, X0 = X0, Y0 = Y0, Z0 = Z0, P = P, S = S)
}

# ------------------------------------------------------------------------------
# 3. Compute G_in covariance matrix for vertex i
# ------------------------------------------------------------------------------
compute_G_in <- function(i, X0, Y0, Z0, tau = 0.05) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  # Compute s_ij = x_0i^T y_0j for all j
  s <- as.vector(X0[i, ] %*% t(Y0))  # n-vector

  # Compute psi'(s_ij)
  w <- dpsi(s, tau = tau)

  # Network contribution: sum_j psi'(s_ij) * y_0j * y_0j^T
  G_net <- matrix(0, d, d)
  for (j in 1:n) {
    G_net <- G_net + w[j] * outer(Y0[j, ], Y0[j, ])
  }

  # Covariate contribution: sum_l z_0l * z_0l^T
  G_cov <- crossprod(Z0)  # Z0^T Z0

  # G_in = (1/(n+p)) * (G_net + G_cov)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

# ------------------------------------------------------------------------------
# 4. Test coverage for a single replicate
# ------------------------------------------------------------------------------
test_single_replicate <- function(n = 50, p_cov = 5, d = 3,
                                   p_sig = 2, q_sig = 1,
                                   vertex_id = 10, alpha = 0.05) {
  cat("=== Testing Single Replicate ===\n")
  cat(sprintf("n = %d, p_cov = %d, d = %d, vertex = %d\n\n",
              n, p_cov, d, vertex_id))

  # Generate data
  cat("1. Generating data...\n")
  data <- generate_grdpg_data(n, p_cov, d, p_sig, q_sig)
  X0 <- data$X0
  Y0 <- data$Y0
  Z0 <- data$Z0

  cat(sprintf("   True position X0[%d,] = (%.3f, %.3f, %.3f)\n",
              vertex_id, X0[vertex_id, 1], X0[vertex_id, 2], X0[vertex_id, 3]))

  # Fit model
  cat("\n2. Fitting GRDPG model...\n")
  fit <- fit_grdpg_cov(data$A, data$B, d = d, p = p_sig, q = q_sig,
                       maxit = 10, tol = 1e-3)

  cat(sprintf("   Converged: %s, Iterations: %d\n",
              fit$converged, fit$iters))
  cat(sprintf("   Estimated position X[%d,] = (%.3f, %.3f, %.3f)\n",
              vertex_id, fit$X[vertex_id, 1], fit$X[vertex_id, 2], fit$X[vertex_id, 3]))

  # Compute covariance matrix G_in^{-1}
  cat("\n3. Computing covariance matrix G_in^{-1}...\n")
  G_in <- compute_G_in(vertex_id, X0, Y0, Z0, tau = fit$tau)

  # Check if G_in is invertible
  eigvals <- eigen(G_in, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("   G_in eigenvalues: min = %.4e, max = %.4e\n",
              min(eigvals), max(eigvals)))

  if (min(eigvals) < 1e-10) {
    cat("   WARNING: G_in is nearly singular!\n")
    return(list(success = FALSE, reason = "singular_G_in"))
  }

  # Invert to get covariance
  Sigma_i <- solve(G_in)
  std_errors <- sqrt(diag(Sigma_i))

  cat(sprintf("   Standard errors: (%.4f, %.4f, %.4f)\n",
              std_errors[1], std_errors[2], std_errors[3]))

  # Note: We need to align estimated X with true X0 via orthogonal rotation
  # For a quick test, let's use Procrustes alignment
  cat("\n4. Aligning estimates with true positions (Procrustes)...\n")

  # Procrustes: find optimal rotation matrix R
  # minimize ||X0 - X_hat R||_F
  M <- t(fit$X) %*% X0
  svd_res <- svd(M)
  R <- svd_res$u %*% t(svd_res$v)

  # Rotate estimated positions
  X_aligned <- fit$X %*% R

  cat(sprintf("   Aligned position X_aligned[%d,] = (%.3f, %.3f, %.3f)\n",
              vertex_id, X_aligned[vertex_id, 1],
              X_aligned[vertex_id, 2], X_aligned[vertex_id, 3]))

  # Compute errors (after rotation)
  error <- X_aligned[vertex_id, ] - X0[vertex_id, ]
  cat(sprintf("   Estimation error: (%.4f, %.4f, %.4f)\n",
              error[1], error[2], error[3]))

  # Check coverage for each coordinate
  cat("\n5. Checking confidence interval coverage...\n")
  z_crit <- qnorm(1 - alpha / 2)  # Two-sided critical value

  coverage <- rep(FALSE, d)
  for (k in 1:d) {
    # CI: x_hat[k] ± z_crit * SE[k]
    # Since we rotated, we need to be careful about the covariance
    # For simplicity in this quick test, use the original SE
    ci_lower <- X_aligned[vertex_id, k] - z_crit * std_errors[k]
    ci_upper <- X_aligned[vertex_id, k] + z_crit * std_errors[k]
    coverage[k] <- (X0[vertex_id, k] >= ci_lower) && (X0[vertex_id, k] <= ci_upper)

    cat(sprintf("   Coordinate %d: CI = [%.4f, %.4f], True = %.4f, Covered = %s\n",
                k, ci_lower, ci_upper, X0[vertex_id, k], coverage[k]))
  }

  cat(sprintf("\n   Overall coverage: %d/%d coordinates\n", sum(coverage), d))

  list(
    success = TRUE,
    coverage = coverage,
    error = error,
    std_errors = std_errors,
    G_in = G_in,
    Sigma_i = Sigma_i,
    fit = fit,
    X0 = X0,
    X_aligned = X_aligned
  )
}

# ------------------------------------------------------------------------------
# RUN THE QUICK TEST
# ------------------------------------------------------------------------------
set.seed(123)

result <- test_single_replicate(
  n = 50,
  p_cov = 5,
  d = 3,
  p_sig = 2,
  q_sig = 1,
  vertex_id = 10,
  alpha = 0.05
)

cat("\n=== Quick Test Complete ===\n")
if (result$success) {
  cat(sprintf("Coverage rate: %d/%d (%.1f%%)\n",
              sum(result$coverage),
              length(result$coverage),
              100 * mean(result$coverage)))
  cat("\nSetup appears viable for full simulation!\n")
} else {
  cat(sprintf("Test failed: %s\n", result$reason))
}
