#!/usr/bin/env Rscript
# Benchmark: Compare original vs optimized G_in computation
# Tests both implementations on same data to measure speedup

library(cgrdpg)

cat("============================================================================\n")
cat("  BENCHMARKING: Original vs Optimized G_in Computation\n")
cat("  Testing on n=200 (representative sample)\n")
cat("============================================================================\n\n")

# Setup
n <- 200
p_cov <- 100
d <- 2
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

cat(sprintf("Test setup: n=%d, p_cov=%d\n\n", n, p_cov))

# ORIGINAL implementation (with for loop)
compute_G_in_original <- function(i, X0, Y0, Z0, tau) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)
  w[i] <- 0

  G_net <- matrix(0, d, d)
  for (j in 1:n) {  # This loop is the bottleneck
    G_net <- G_net + w[j] * outer(Y0[j, ], Y0[j, ])
  }

  G_cov <- crossprod(Z0)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

# OPTIMIZED implementation (vectorized)
compute_G_in_optimized <- function(i, X0, Y0, Z0, tau) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)
  w[i] <- 0

  # Vectorized: t(Y) %*% diag(w) %*% Y
  Y_weighted <- Y0 * sqrt(w)
  G_net <- crossprod(Y_weighted)

  G_cov <- crossprod(Z0)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

# Test correctness first
cat("Testing correctness (comparing outputs)...\n")
test_vertex <- 50
G_orig <- compute_G_in_original(test_vertex, X0, Y0, Z0, tau)
G_opt <- compute_G_in_optimized(test_vertex, X0, Y0, Z0, tau)
max_diff <- max(abs(G_orig - G_opt))
cat(sprintf("  Max difference between methods: %.2e\n", max_diff))
if (max_diff < 1e-10) {
  cat("  ✓ Methods produce identical results!\n\n")
} else {
  cat("  ✗ WARNING: Methods differ!\n\n")
}

# Benchmark on multiple vertices
n_vertices_to_test <- 50
vertices_to_test <- sample(1:n, n_vertices_to_test)

cat(sprintf("Benchmarking on %d random vertices...\n\n", n_vertices_to_test))

# Warm up
for (i in 1:5) {
  G_orig <- compute_G_in_original(i, X0, Y0, Z0, tau)
  G_opt <- compute_G_in_optimized(i, X0, Y0, Z0, tau)
}

# Benchmark ORIGINAL
cat("Running ORIGINAL implementation...\n")
time_orig_start <- Sys.time()
for (i in vertices_to_test) {
  G_orig <- compute_G_in_original(i, X0, Y0, Z0, tau)
}
time_orig <- as.numeric(difftime(Sys.time(), time_orig_start, units = "secs"))
cat(sprintf("  Time: %.3f seconds (%.1f ms per vertex)\n\n",
            time_orig, 1000 * time_orig / n_vertices_to_test))

# Benchmark OPTIMIZED
cat("Running OPTIMIZED implementation...\n")
time_opt_start <- Sys.time()
for (i in vertices_to_test) {
  G_opt <- compute_G_in_optimized(i, X0, Y0, Z0, tau)
}
time_opt <- as.numeric(difftime(Sys.time(), time_opt_start, units = "secs"))
cat(sprintf("  Time: %.3f seconds (%.1f ms per vertex)\n\n",
            time_opt, 1000 * time_opt / n_vertices_to_test))

# Results
speedup <- time_orig / time_opt
cat("============================================================================\n")
cat("  BENCHMARK RESULTS\n")
cat("============================================================================\n\n")
cat(sprintf("Original:  %.3f seconds\n", time_orig))
cat(sprintf("Optimized: %.3f seconds\n", time_opt))
cat(sprintf("Speedup:   %.2fx faster\n\n", speedup))

# Extrapolate to full simulation
cat("Estimated time for full n=1000 simulation (all vertices):\n")
time_per_vertex_orig <- time_orig / n_vertices_to_test
time_per_vertex_opt <- time_opt / n_vertices_to_test
n_full <- 1000

est_time_orig <- time_per_vertex_orig * n_full * 2  # x2 for TRUE + PLUGIN
est_time_opt <- time_per_vertex_opt * n_full * 2

cat(sprintf("  Original:  %.1f minutes\n", est_time_orig / 60))
cat(sprintf("  Optimized: %.1f minutes\n", est_time_opt / 60))
cat(sprintf("  Time saved: %.1f minutes per replication\n", (est_time_orig - est_time_opt) / 60))
cat(sprintf("  For 100 reps: %.1f hours saved!\n\n", (est_time_orig - est_time_opt) * 100 / 3600))

# Save benchmark results
benchmark_results <- list(
  n = n,
  p_cov = p_cov,
  n_vertices_tested = n_vertices_to_test,
  time_original = time_orig,
  time_optimized = time_opt,
  speedup = speedup,
  max_difference = max_diff,
  estimated_time_n1000_original_mins = est_time_orig / 60,
  estimated_time_n1000_optimized_mins = est_time_opt / 60
)

saveRDS(benchmark_results, "benchmark_results.rds")
cat("Benchmark results saved to: benchmark_results.rds\n")

cat("============================================================================\n")
cat("\nRECOMMENDATION:\n")
if (speedup > 3) {
  cat(sprintf("The optimized version is %.1fx faster! Highly recommended to use.\n", speedup))
} else if (speedup > 1.5) {
  cat(sprintf("The optimized version is %.1fx faster. Recommended for large simulations.\n", speedup))
} else {
  cat("The speedup is modest. May not be worth switching.\n")
}
cat("============================================================================\n")
