#!/usr/bin/env Rscript
# Profile the vertex-wise coverage computation to identify bottlenecks
# Uses a smaller n to make profiling faster while still being representative

library(cgrdpg)

cat("============================================================================\n")
cat("  PROFILING VERTEX-WISE COVERAGE COMPUTATION\n")
cat("  Using reduced n=200 for faster profiling\n")
cat("============================================================================\n\n")

# Reduced parameters for faster profiling
n <- 200
p_cov <- 100
d <- 2
maxit <- 30
tol <- 0.01
tau <- 0.005
set.seed(598)

cat("Parameters: n=%d, p_cov=%d, d=%d, tau=%.3f\n\n", n, p_cov, d, tau)

# Generate data (same setup as main script)
cat("Generating data...\n")
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

cat("Fitting model...\n")
fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                     maxit = maxit, tol = tol, tau = tau)

# Oracle Procrustes alignment
X_est <- fit$X
M <- t(X_est) %*% X0
svd_res <- svd(M)
Q <- svd_res$u %*% t(svd_res$v)
X_aligned <- X_est %*% Q

Z_est <- fit$Z
Y_aligned <- X_aligned %*% S
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

# Define the G_in computation functions
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

cat("\n============================================================================\n")
cat("  STARTING PROFILING (this will take a few moments...)\n")
cat("============================================================================\n\n")

# Profile the coverage computation
chi2_crit <- qchisq(0.95, df = d)
in_ellipse_true <- logical(n)
in_ellipse_plugin <- logical(n)

# Start profiling
Rprof("profile_output.out", memory.profiling = TRUE, interval = 0.01)

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

Rprof(NULL)

cat("Profiling complete!\n\n")

# Analyze profiling results
cat("============================================================================\n")
cat("  PROFILING RESULTS\n")
cat("============================================================================\n\n")

prof_summary <- summaryRprof("profile_output.out", memory = "both")

cat("--- TIME PROFILING (by function) ---\n")
print(head(prof_summary$by.self, 15))

cat("\n--- TIME PROFILING (by total) ---\n")
print(head(prof_summary$by.total, 15))

cat("\n--- MEMORY PROFILING ---\n")
if (!is.null(prof_summary$by.self.memory)) {
  print(head(prof_summary$by.self.memory, 10))
}

# Save detailed results
save(prof_summary, file = "profiling_results.rds")
cat("\n\nDetailed results saved to: profiling_results.rds\n")

# Cleanup
file.remove("profile_output.out")

cat("\n============================================================================\n")
cat("  KEY BOTTLENECKS IDENTIFIED\n")
cat("============================================================================\n\n")

cat("The profiling output above shows:\n")
cat("  - 'by.self': Time spent IN that function (excluding calls to other functions)\n")
cat("  - 'by.total': Total time including all nested function calls\n\n")

cat("Look for:\n")
cat("  1. Functions with high 'self.time' - these are direct bottlenecks\n")
cat("  2. Functions called many times - candidates for optimization\n")
cat("  3. Matrix operations (outer, crossprod, %*%, solve) - potential for vectorization\n")
cat("  4. The dpsi() function - may be called O(n^2) times\n\n")

cat("Recommendations based on typical bottlenecks:\n")
cat("  - If dpsi() is slow: Consider vectorizing or caching results\n")
cat("  - If outer() loops are slow: Consider pre-computing or using matrix operations\n")
cat("  - If eigen() is slow: Consider using only.values=TRUE (already done)\n")
cat("  - If solve() is slow: Consider caching matrix inverse\n\n")

cat("============================================================================\n")
