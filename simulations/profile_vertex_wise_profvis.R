#!/usr/bin/env Rscript
# Profile using profvis for interactive visualization
# Install profvis with: install.packages("profvis")

if (!requireNamespace("profvis", quietly = TRUE)) {
  cat("Installing profvis package...\n")
  install.packages("profvis", repos = "https://cloud.r-project.org")
}

library(cgrdpg)
library(profvis)

cat("============================================================================\n")
cat("  PROFVIS PROFILING - Vertex-wise Coverage\n")
cat("  Using n=200 for representative profiling\n")
cat("============================================================================\n\n")

# Profiling with profvis
p <- profvis({
  # Setup
  n <- 200
  p_cov <- 100
  d <- 2
  maxit <- 30
  tol <- 0.01
  tau <- 0.005
  set.seed(598)

  # Generate data
  i_vals <- 1:n
  theta <- pi * i_vals / (n - 1)
  X0 <- matrix(0, n, 2)
  X0[, 1] <- 0.28 * sin(theta) + 0.42
  X0[, 2] <- 0.28 * cos(theta) + 0.42

  S <- diag(c(1, 1))
  Y0 <- X0 %*% S
  Z0 <- matrix(rnorm(p_cov * d), p_cov, d)
  P <- X0 %*% t(Y0)

  A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
  A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

  # Fit model
  fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                       maxit = maxit, tol = tol, tau = tau)

  # Alignment
  X_est <- fit$X
  M <- t(X_est) %*% X0
  svd_res <- svd(M)
  Q <- svd_res$u %*% t(svd_res$v)
  X_aligned <- X_est %*% Q

  Z_est <- fit$Z
  Y_aligned <- X_aligned %*% S
  Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

  # G_in computation functions
  compute_G_in_true <- function(i, X0, Y0, Z0, tau) {
    n <- nrow(X0)
    p_cov <- nrow(Z0)
    d <- ncol(X0)

    s <- as.vector(X0[i, ] %*% t(Y0))
    w <- dpsi(s, tau = tau)
    w[i] <- 0

    G_net <- matrix(0, d, d)
    for (j in 1:n) {
      G_net <- G_net + w[j] * outer(Y0[j, ], Y0[j, ])
    }

    G_cov <- crossprod(Z0)
    G_in <- (G_net + G_cov) / (n + p_cov)

    return(G_in)
  }

  compute_G_in_plugin <- function(i, X_est, Y_est, Z_est, tau) {
    n <- nrow(X_est)
    p_cov <- nrow(Z_est)
    d <- ncol(X_est)

    s <- as.vector(X_est[i, ] %*% t(Y_est))
    w <- dpsi(s, tau = tau)
    w[i] <- 0

    G_net <- matrix(0, d, d)
    for (j in 1:n) {
      G_net <- G_net + w[j] * outer(Y_est[j, ], Y_est[j, ])
    }

    G_cov <- crossprod(Z_est)
    G_in <- (G_net + G_cov) / (n + p_cov)

    return(G_in)
  }

  # Coverage computation
  chi2_crit <- qchisq(0.95, df = d)
  in_ellipse_true <- logical(n)
  in_ellipse_plugin <- logical(n)

  for (i in 1:n) {
    # TRUE G_in
    G_in_true <- compute_G_in_true(i, X0, Y0, Z0, tau = tau)
    eig_vals <- eigen(G_in_true, only.values = TRUE)$values
    if (min(eig_vals) >= 1e-10) {
      diff_vec <- X0[i, ] - X_aligned[i, ]
      mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
      in_ellipse_true[i] <- (mahal_dist <= chi2_crit)
    } else {
      in_ellipse_true[i] <- NA
    }

    # PLUG-IN G_in
    G_in_plugin <- compute_G_in_plugin(i, X_aligned, Y_aligned, Z_updated, tau = tau)
    eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
    if (min(eig_vals) >= 1e-10) {
      diff_vec <- X0[i, ] - X_aligned[i, ]
      mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_plugin %*% diff_vec)
      in_ellipse_plugin[i] <- (mahal_dist <= chi2_crit)
    } else {
      in_ellipse_plugin[i] <- NA
    }
  }
})

# Save the profvis output
htmlwidgets::saveWidget(p, "profiling_visualization.html")

cat("\n============================================================================\n")
cat("  PROFVIS OUTPUT SAVED\n")
cat("============================================================================\n\n")
cat("Interactive profiling visualization saved to: profiling_visualization.html\n")
cat("Open this file in a web browser to see:\n")
cat("  - Flame graph showing function call hierarchy\n")
cat("  - Time spent in each function\n")
cat("  - Line-by-line profiling of your code\n\n")

cat("To view the results:\n")
cat("  open profiling_visualization.html\n\n")

# Also print text summary
cat("Text summary of profiling:\n")
print(p)
