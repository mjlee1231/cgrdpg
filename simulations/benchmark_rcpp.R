#!/usr/bin/env Rscript
# Benchmark: R vs Rcpp implementations
# Tests speedup from C++ versions of psi functions and Fisher scoring

library(cgrdpg)

cat("============================================================================\n")
cat("  BENCHMARK: R vs Rcpp Performance\n")
cat("  Testing C++ optimizations for key bottlenecks\n")
cat("============================================================================\n\n")

# Check if Rcpp functions are available
rcpp_available <- exists("dpsi_cpp") && exists("compute_G_in_vectorized_cpp")

if (!rcpp_available) {
  cat("⚠ Rcpp functions not found! Package needs to be compiled with Rcpp support.\n")
  cat("Run: R CMD INSTALL . (with Rcpp installed)\n\n")
  quit(status = 1)
}

# Setup test data
n <- 500
p_cov <- 250
d <- 2
tau <- 0.005
set.seed(598)

cat(sprintf("Test setup: n=%d, p_cov=%d, d=%d\n\n", n, p_cov, d))

# Generate latent positions
theta <- pi * (1:n) / (n - 1)
X0 <- matrix(0, n, d)
X0[, 1] <- 0.28 * sin(theta) + 0.42
X0[, 2] <- 0.28 * cos(theta) + 0.42
Y0 <- X0
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

# Test vector
s <- as.vector(X0 %*% t(Y0))[1:n]

cat("============================================================================\n")
cat("  TEST 1: psi Functions\n")
cat("============================================================================\n\n")

# Benchmark psi
n_reps <- 1000
cat(sprintf("Running psi() %d times...\n", n_reps))

time_r_start <- Sys.time()
for (i in 1:n_reps) {
  result_r <- psi(s, tau)
}
time_r_psi <- as.numeric(difftime(Sys.time(), time_r_start, units = "secs"))

time_cpp_start <- Sys.time()
for (i in 1:n_reps) {
  result_cpp <- psi_cpp(s, tau)
}
time_cpp_psi <- as.numeric(difftime(Sys.time(), time_cpp_start, units = "secs"))

# Check correctness
max_diff_psi <- max(abs(result_r - result_cpp))

cat(sprintf("  R version:   %.4f sec\n", time_r_psi))
cat(sprintf("  C++ version: %.4f sec\n", time_cpp_psi))
cat(sprintf("  Speedup:     %.2fx\n", time_r_psi / time_cpp_psi))
cat(sprintf("  Max diff:    %.2e (should be ~0)\n\n", max_diff_psi))

# Benchmark dpsi
cat(sprintf("Running dpsi() %d times...\n", n_reps))

time_r_start <- Sys.time()
for (i in 1:n_reps) {
  result_r <- dpsi(s, tau)
}
time_r_dpsi <- as.numeric(difftime(Sys.time(), time_r_start, units = "secs"))

time_cpp_start <- Sys.time()
for (i in 1:n_reps) {
  result_cpp <- dpsi_cpp(s, tau)
}
time_cpp_dpsi <- as.numeric(difftime(Sys.time(), time_cpp_start, units = "secs"))

max_diff_dpsi <- max(abs(result_r - result_cpp))

cat(sprintf("  R version:   %.4f sec\n", time_r_dpsi))
cat(sprintf("  C++ version: %.4f sec\n", time_cpp_dpsi))
cat(sprintf("  Speedup:     %.2fx\n", time_r_dpsi / time_cpp_dpsi))
cat(sprintf("  Max diff:    %.2e (should be ~0)\n\n", max_diff_dpsi))

# Benchmark Psi
cat(sprintf("Running Psi() %d times...\n", n_reps))

time_r_start <- Sys.time()
for (i in 1:n_reps) {
  result_r <- Psi(s, tau)
}
time_r_Psi <- as.numeric(difftime(Sys.time(), time_r_start, units = "secs"))

time_cpp_start <- Sys.time()
for (i in 1:n_reps) {
  result_cpp <- Psi_cpp(s, tau)
}
time_cpp_Psi <- as.numeric(difftime(Sys.time(), time_cpp_start, units = "secs"))

max_diff_Psi <- max(abs(result_r - result_cpp))

cat(sprintf("  R version:   %.4f sec\n", time_r_Psi))
cat(sprintf("  C++ version: %.4f sec\n", time_cpp_Psi))
cat(sprintf("  Speedup:     %.2fx\n", time_r_Psi / time_cpp_Psi))
cat(sprintf("  Max diff:    %.2e (should be ~0)\n\n", max_diff_Psi))

cat("============================================================================\n")
cat("  TEST 2: G_in Computation (Coverage Bottleneck)\n")
cat("============================================================================\n\n")

# R version (vectorized)
compute_G_in_r <- function(i, X, Y, Z, tau) {
  n <- nrow(X)
  p_cov <- nrow(Z)

  s <- as.vector(X[i, ] %*% t(Y))
  w <- dpsi(s, tau = tau)
  w[i] <- 0

  Y_weighted <- Y * sqrt(w)
  G_net <- crossprod(Y_weighted)
  G_cov <- crossprod(Z)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

n_vertices <- 50
cat(sprintf("Computing G_in for %d vertices...\n", n_vertices))

time_r_start <- Sys.time()
for (i in 1:n_vertices) {
  G_r <- compute_G_in_r(i, X0, Y0, Z0, tau)
}
time_r_G <- as.numeric(difftime(Sys.time(), time_r_start, units = "secs"))

time_cpp_start <- Sys.time()
for (i in 1:n_vertices) {
  G_cpp <- compute_G_in_vectorized_cpp(i - 1, X0, Y0, Z0, tau)  # C++ is 0-indexed
}
time_cpp_G <- as.numeric(difftime(Sys.time(), time_cpp_start, units = "secs"))

max_diff_G <- max(abs(G_r - G_cpp))

cat(sprintf("  R version (vectorized): %.4f sec\n", time_r_G))
cat(sprintf("  C++ version:            %.4f sec\n", time_cpp_G))
cat(sprintf("  Speedup:                %.2fx\n", time_r_G / time_cpp_G))
cat(sprintf("  Max diff:               %.2e (should be ~0)\n\n", max_diff_G))

cat("============================================================================\n")
cat("  SUMMARY\n")
cat("============================================================================\n\n")

avg_speedup_psi <- mean(c(
  time_r_psi / time_cpp_psi,
  time_r_dpsi / time_cpp_dpsi,
  time_r_Psi / time_cpp_Psi
))

cat(sprintf("Average speedup (psi functions): %.2fx\n", avg_speedup_psi))
cat(sprintf("Speedup (G_in computation):      %.2fx\n", time_r_G / time_cpp_G))
cat(sprintf("\nAll correctness checks passed: %s\n",
            ifelse(max(max_diff_psi, max_diff_dpsi, max_diff_Psi, max_diff_G) < 1e-10,
                   "✓ YES", "✗ NO")))

cat("\n============================================================================\n")
cat("  ESTIMATED IMPACT ON FULL SIMULATION\n")
cat("============================================================================\n\n")

cat("For n=2000 simulation:\n")
cat("  - dpsi() called ~4M times per replication\n")
cat(sprintf("  - With %.2fx speedup: Saves ~5-10%% of total time\n", avg_speedup_psi))
cat(sprintf("  - G_in speedup: %.2fx -> Saves ~10-20%% of coverage time\n", time_r_G / time_cpp_G))
cat("  - Combined with parallel: Total ~4-6x speedup over serial\n\n")

cat("Recommendation:\n")
if (avg_speedup_psi > 5 && time_r_G / time_cpp_G > 2) {
  cat("  ✓ HIGHLY RECOMMENDED: Significant Rcpp speedup achieved!\n")
  cat("    Deploy to HPC for maximum performance.\n")
} else if (avg_speedup_psi > 2) {
  cat("  ✓ RECOMMENDED: Modest speedup, worthwhile for large simulations.\n")
} else {
  cat("  ⚠ Limited speedup - parallel optimization may be sufficient.\n")
}

cat("\n============================================================================\n")

# Save results
benchmark_results <- list(
  speedup_psi = time_r_psi / time_cpp_psi,
  speedup_dpsi = time_r_dpsi / time_cpp_dpsi,
  speedup_Psi = time_r_Psi / time_cpp_Psi,
  speedup_G_in = time_r_G / time_cpp_G,
  avg_speedup_psi_functions = avg_speedup_psi,
  max_diff_psi = max_diff_psi,
  max_diff_dpsi = max_diff_dpsi,
  max_diff_Psi = max_diff_Psi,
  max_diff_G = max_diff_G
)

saveRDS(benchmark_results, "benchmark_rcpp_results.rds")
cat("Benchmark results saved to: benchmark_rcpp_results.rds\n")
