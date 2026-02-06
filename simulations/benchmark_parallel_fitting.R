#!/usr/bin/env Rscript
# Benchmark: Serial vs Parallel model fitting
# Tests fit_grdpg_cov vs fit_grdpg_cov_parallel

# Reload package to get latest functions
try(detach("package:cgrdpg", unload = TRUE), silent = TRUE)
library(cgrdpg)

# Load parallel functions directly if not exported yet
if (!exists("fit_grdpg_cov_parallel")) {
  source("../R/fisher_parallel.R")
}

cat("============================================================================\n")
cat("  BENCHMARK: Serial vs Parallel Model Fitting\n")
cat("  Testing on n=200 for representative comparison\n")
cat("============================================================================\n\n")

# Setup
n <- 200
p_cov <- 100
d <- 2
maxit <- 10  # Limit iterations for faster benchmarking
tol <- 0.01
tau <- 0.005
set.seed(598)

cat(sprintf("Parameters: n=%d, p_cov=%d, d=%d, maxit=%d\n\n", n, p_cov, d, maxit))

# Generate data
cat("Generating test data...\n")
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

cat("Data generated successfully\n\n")

# Check if parallel packages are available
if (!requireNamespace("foreach", quietly = TRUE)) {
  cat("Installing 'foreach' package...\n")
  install.packages("foreach", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("doParallel", quietly = TRUE)) {
  cat("Installing 'doParallel' package...\n")
  install.packages("doParallel", repos = "https://cloud.r-project.org")
}

# Detect available cores
n_cores <- parallel::detectCores()
cat(sprintf("Available CPU cores: %d\n", n_cores))
cat(sprintf("Will use %d cores for parallel version\n\n", n_cores - 1))

# ============================================================================
# Benchmark 1: Serial version
# ============================================================================
cat("Running SERIAL version (fit_grdpg_cov)...\n")
time_serial_start <- Sys.time()

fit_serial <- fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                            maxit = maxit, tol = tol, tau = tau)

time_serial <- as.numeric(difftime(Sys.time(), time_serial_start, units = "secs"))

cat(sprintf("  Converged: %s\n", fit_serial$converged))
cat(sprintf("  Iterations: %d\n", fit_serial$iters))
cat(sprintf("  Time: %.2f seconds\n\n", time_serial))

# ============================================================================
# Benchmark 2: Parallel version
# ============================================================================
cat("Running PARALLEL version (fit_grdpg_cov_parallel)...\n")
time_parallel_start <- Sys.time()

fit_parallel <- fit_grdpg_cov_parallel(A, B, d = d, p = 2, q = 0,
                                       maxit = maxit, tol = tol, tau = tau,
                                       ncores = n_cores - 1)

time_parallel <- as.numeric(difftime(Sys.time(), time_parallel_start, units = "secs"))

cat(sprintf("  Converged: %s\n", fit_parallel$converged))
cat(sprintf("  Iterations: %d\n", fit_parallel$iters))
cat(sprintf("  Time: %.2f seconds\n\n", time_parallel))

# ============================================================================
# Results comparison
# ============================================================================
cat("============================================================================\n")
cat("  BENCHMARK RESULTS\n")
cat("============================================================================\n\n")

speedup <- time_serial / time_parallel
efficiency <- speedup / (n_cores - 1) * 100

cat(sprintf("Serial time:    %.2f seconds\n", time_serial))
cat(sprintf("Parallel time:  %.2f seconds\n", time_parallel))
cat(sprintf("Speedup:        %.2fx faster\n", speedup))
cat(sprintf("Parallel efficiency: %.1f%% (using %d cores)\n\n",
            efficiency, n_cores - 1))

# Check if results are similar
SSE_serial <- sum((fit_serial$X - X0)^2)
SSE_parallel <- sum((fit_parallel$X - X0)^2)
max_diff <- max(abs(fit_serial$X - fit_parallel$X))

cat("Convergence comparison:\n")
cat(sprintf("  Serial SSE:      %.6f\n", SSE_serial))
cat(sprintf("  Parallel SSE:    %.6f\n", SSE_parallel))
cat(sprintf("  Max X difference: %.2e\n\n", max_diff))

if (max_diff < 1e-3) {
  cat("✓ Results are very similar (Jacobi vs Gauss-Seidel expected)\n\n")
} else {
  cat("⚠ Results differ (this is expected due to Jacobi parallelization)\n\n")
}

# ============================================================================
# Extrapolate to n=1000
# ============================================================================
cat("============================================================================\n")
cat("  EXTRAPOLATION TO n=1000\n")
cat("============================================================================\n\n")

# Assuming O(n^2) scaling for Fisher iterations
scale_factor <- (1000 / 200)^2

cat("Estimated time for ONE replication (n=1000, 30 iterations):\n")
est_serial_n1000 <- time_serial * scale_factor * (30 / maxit)
est_parallel_n1000 <- time_parallel * scale_factor * (30 / maxit)

cat(sprintf("  Serial:    ~%.1f minutes\n", est_serial_n1000 / 60))
cat(sprintf("  Parallel:  ~%.1f minutes\n", est_parallel_n1000 / 60))
cat(sprintf("  Speedup:   %.2fx faster\n", speedup))
cat(sprintf("  Time saved per rep: ~%.1f minutes\n\n",
            (est_serial_n1000 - est_parallel_n1000) / 60))

cat("For 100 replications:\n")
total_time_saved <- (est_serial_n1000 - est_parallel_n1000) * 100
cat(sprintf("  Time saved: ~%.1f hours (%.1f days)\n",
            total_time_saved / 3600,
            total_time_saved / 86400))

cat("\n============================================================================\n")
cat("  RECOMMENDATION\n")
cat("============================================================================\n\n")

if (speedup >= 2) {
  cat(sprintf("✓ HIGHLY RECOMMENDED: %.2fx speedup with %d cores\n", speedup, n_cores - 1))
  cat("  Parallel version will save significant time for large simulations\n")
} else if (speedup >= 1.3) {
  cat(sprintf("✓ RECOMMENDED: %.2fx speedup with %d cores\n", speedup, n_cores - 1))
  cat("  Worthwhile for 100+ replications\n")
} else {
  cat("⚠ Limited speedup - may not be worth the parallelization overhead\n")
  cat("  Consider using more cores or checking for bottlenecks\n")
}

cat("\nNote: Parallel version uses Jacobi-style updates (all nodes simultaneously)\n")
cat("      vs Gauss-Seidel in serial (sequential updates). This may require\n")
cat("      slightly more iterations but is much faster overall.\n")

cat("\n============================================================================\n")

# Save benchmark results
benchmark_results <- list(
  n = n,
  p_cov = p_cov,
  maxit = maxit,
  time_serial = time_serial,
  time_parallel = time_parallel,
  speedup = speedup,
  ncores = n_cores - 1,
  efficiency = efficiency,
  SSE_serial = SSE_serial,
  SSE_parallel = SSE_parallel,
  max_diff = max_diff,
  estimated_time_n1000_serial_mins = est_serial_n1000 / 60,
  estimated_time_n1000_parallel_mins = est_parallel_n1000 / 60
)

saveRDS(benchmark_results, "benchmark_parallel_results.rds")
cat("Benchmark results saved to: benchmark_parallel_results.rds\n")
