#!/usr/bin/env Rscript
# Comprehensive benchmark: Full simulation workflow timing breakdown
# Identifies what proportion of time is spent in each step

library(cgrdpg)

cat("============================================================================\n")
cat("  FULL WORKFLOW TIMING ANALYSIS\n")
cat("  Identifying all computational bottlenecks (n=200)\n")
cat("============================================================================\n\n")

n <- 200
p_cov <- 100
d <- 2
maxit <- 30
tol <- 0.01
tau <- 0.005
set.seed(598)

# Timing storage
timings <- list()

# ============================================================================
# 1. Data Generation
# ============================================================================
cat("1. Generating data...\n")
t_start <- Sys.time()

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

timings$data_generation <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("   Time: %.3f seconds\n\n", timings$data_generation))

# ============================================================================
# 2. Model Fitting
# ============================================================================
cat("2. Fitting GRDPG model...\n")
t_start <- Sys.time()

fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                     maxit = maxit, tol = tol, tau = tau)

timings$model_fitting <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("   Time: %.3f seconds\n\n", timings$model_fitting))

# ============================================================================
# 3. Procrustes Alignment
# ============================================================================
cat("3. Procrustes alignment...\n")
t_start <- Sys.time()

X_est <- fit$X
M <- t(X_est) %*% X0
svd_res <- svd(M)
Q <- svd_res$u %*% t(svd_res$v)
X_aligned <- X_est %*% Q

Z_est <- fit$Z
Y_aligned <- X_aligned %*% S
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

timings$alignment <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("   Time: %.3f seconds\n\n", timings$alignment))

# ============================================================================
# 4. Coverage Computation (ORIGINAL - with loop)
# ============================================================================
cat("4. Coverage computation (ORIGINAL with for loop)...\n")

compute_G_in_original <- function(i, X0, Y0, Z0, tau) {
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

chi2_crit <- qchisq(0.95, df = d)
in_ellipse_true <- logical(n)
in_ellipse_plugin <- logical(n)

t_start <- Sys.time()

for (i in 1:n) {
  # TRUE G_in
  G_in_true <- compute_G_in_original(i, X0, Y0, Z0, tau = tau)
  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
    in_ellipse_true[i] <- (mahal_dist <= chi2_crit)
  }

  # PLUG-IN G_in
  G_in_plugin <- compute_G_in_original(i, X_aligned, Y_aligned, Z_updated, tau = tau)
  eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_plugin %*% diff_vec)
    in_ellipse_plugin[i] <- (mahal_dist <= chi2_crit)
  }
}

timings$coverage_original <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("   Time: %.3f seconds\n\n", timings$coverage_original))

# ============================================================================
# 5. Coverage Computation (OPTIMIZED - vectorized)
# ============================================================================
cat("5. Coverage computation (OPTIMIZED vectorized)...\n")

compute_G_in_optimized <- function(i, X0, Y0, Z0, tau) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)
  w[i] <- 0

  Y_weighted <- Y0 * sqrt(w)
  G_net <- crossprod(Y_weighted)

  G_cov <- crossprod(Z0)
  G_in <- (G_net + G_cov) / (n + p_cov)
  return(G_in)
}

in_ellipse_true_opt <- logical(n)
in_ellipse_plugin_opt <- logical(n)

t_start <- Sys.time()

for (i in 1:n) {
  # TRUE G_in
  G_in_true <- compute_G_in_optimized(i, X0, Y0, Z0, tau = tau)
  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
    in_ellipse_true_opt[i] <- (mahal_dist <= chi2_crit)
  }

  # PLUG-IN G_in
  G_in_plugin <- compute_G_in_optimized(i, X_aligned, Y_aligned, Z_updated, tau = tau)
  eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_plugin %*% diff_vec)
    in_ellipse_plugin_opt[i] <- (mahal_dist <= chi2_crit)
  }
}

timings$coverage_optimized <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("   Time: %.3f seconds\n\n", timings$coverage_optimized))

# ============================================================================
# Summary
# ============================================================================
cat("============================================================================\n")
cat("  TIMING BREAKDOWN (n=200)\n")
cat("============================================================================\n\n")

total_original <- sum(unlist(timings[c("data_generation", "model_fitting", "alignment", "coverage_original")]))
total_optimized <- sum(unlist(timings[c("data_generation", "model_fitting", "alignment", "coverage_optimized")]))

cat(sprintf("%-25s %8.3f sec (%5.1f%% of total)\n",
            "1. Data Generation:",
            timings$data_generation,
            100 * timings$data_generation / total_original))

cat(sprintf("%-25s %8.3f sec (%5.1f%% of total) <- MAIN BOTTLENECK\n",
            "2. Model Fitting:",
            timings$model_fitting,
            100 * timings$model_fitting / total_original))

cat(sprintf("%-25s %8.3f sec (%5.1f%% of total)\n",
            "3. Procrustes Alignment:",
            timings$alignment,
            100 * timings$alignment / total_original))

cat(sprintf("%-25s %8.3f sec (%5.1f%% of total)\n",
            "4. Coverage (ORIGINAL):",
            timings$coverage_original,
            100 * timings$coverage_original / total_original))

cat(sprintf("%-25s %8.3f sec (%5.1f%% of total)\n",
            "5. Coverage (OPTIMIZED):",
            timings$coverage_optimized,
            100 * timings$coverage_optimized / total_optimized))

cat(sprintf("\n%-25s %8.3f sec\n", "TOTAL (Original):", total_original))
cat(sprintf("%-25s %8.3f sec\n", "TOTAL (Optimized):", total_optimized))
cat(sprintf("%-25s %8.3f sec (%.1f%% faster overall)\n\n",
            "Time Saved:",
            total_original - total_optimized,
            100 * (total_original - total_optimized) / total_original))

cat("Coverage speedup: %.1fx faster\n\n", timings$coverage_original / timings$coverage_optimized)

# Extrapolate to n=1000
cat("============================================================================\n")
cat("  EXTRAPOLATION TO n=1000 (full simulation)\n")
cat("============================================================================\n\n")

scale_factor <- (1000 / 200)^2  # Many operations scale as O(n^2)

cat("Estimated times for ONE replication (n=1000):\n")
cat(sprintf("  Model fitting:         ~%.1f minutes\n",
            timings$model_fitting * scale_factor / 60))
cat(sprintf("  Coverage (ORIGINAL):   ~%.1f minutes\n",
            timings$coverage_original * scale_factor / 60))
cat(sprintf("  Coverage (OPTIMIZED):  ~%.1f minutes\n",
            timings$coverage_optimized * scale_factor / 60))
cat(sprintf("  Time saved per rep:    ~%.1f minutes\n\n",
            (timings$coverage_original - timings$coverage_optimized) * scale_factor / 60))

cat("For 100 replications:\n")
cat(sprintf("  Time saved:            ~%.1f hours\n",
            (timings$coverage_original - timings$coverage_optimized) * scale_factor * 100 / 3600))

cat("\n============================================================================\n")
cat("  RECOMMENDATIONS\n")
cat("============================================================================\n\n")

cat("1. MAIN BOTTLENECK: Model fitting (fit_grdpg_cov)\n")
cat(sprintf("   Takes %.1f%% of total time - limited room for optimization here\n\n",
            100 * timings$model_fitting / total_original))

cat("2. Coverage computation speedup: %.1fx with vectorization\n",
    timings$coverage_original / timings$coverage_optimized)
cat("   Recommended to use OPTIMIZED version\n\n")

cat("3. Overall speedup: %.1f%%\n",
    100 * (total_original - total_optimized) / total_original)
cat("   Modest but worthwhile for 100 replications\n\n")

cat("============================================================================\n")

# Save results
saveRDS(timings, "full_workflow_timings.rds")
cat("\nDetailed timings saved to: full_workflow_timings.rds\n")
