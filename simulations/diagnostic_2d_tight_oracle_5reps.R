# Diagnostic: 2D TIGHT latent positions with Manual Oracle Procrustes Alignment
# 5 replications, ALL nodes checked
# NOTE: Package does NOT include Oracle Procrustes - alignment is done manually in this script
# ARRAY JOB VERSION: Accepts replication number from command line
# PARALLEL VERSION: Uses mclapply for parallel coverage computation
library(cgrdpg)
library(parallel)

# Get replication number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript diagnostic_2d_tight_oracle_5reps.R <rep_number>")
}
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 5) {
  stop("Rep number must be between 1 and 5")
}

# Get number of cores from SLURM environment
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores) || n_cores < 1) {
  n_cores <- 1  # fallback for local runs
  cat("Warning: SLURM_CPUS_PER_TASK not set, using 1 core\n")
} else {
  cat(sprintf("Using %d cores for parallel processing\n", n_cores))
}

compute_G_in_true <- function(i, X0, Y0, Z0, tau) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)
  w[i] <- 0  # Exclude diagonal (no self-loops)

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
  w[i] <- 0  # Exclude diagonal (no self-loops)

  G_net <- matrix(0, d, d)
  for (j in 1:n) {
    G_net <- G_net + w[j] * outer(Y_est[j, ], Y_est[j, ])
  }

  G_cov <- crossprod(Z_est)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

cat("============================================================================\n")
cat("  DIAGNOSTIC: 2D TIGHT - Manual Oracle Procrustes Alignment\n")
cat(sprintf("  Running REPLICATION %d (Array Job)\n", rep_id))
cat(sprintf("  ALL 500 nodes checked with %d parallel cores\n", n_cores))
cat("  Package defaults: tau=0.001, ls_beta=0.35, tol=0.005, maxit=30\n")
cat("  Script overrides: tau=0.005, tol=0.01, ADAPTIVE ls_beta=0.8->0.4\n")
cat("  With diagonal augmentation for ASE initialization\n")
cat("============================================================================\n\n")

# Fixed parameters
n <- 500
p_cov <- 250
d <- 2
maxit <- 30
tol <- 0.01

# Set seed based on replication number for reproducibility
set.seed(123 + rep_id)

# Generate 2D latent positions (TIGHT configuration)
# Target Max Prob: ~0.91
cat("Generating 2D TIGHT latent positions...\n")
i_vals <- 0:(n-1)
theta <- pi * i_vals / (n - 1)
X0 <- matrix(0, n, 2)
X0[, 1] <- 0.12 * sin(theta) + 0.55
X0[, 2] <- 0.12 * cos(theta) + 0.55

S <- diag(c(1, 1))
Y0 <- X0 %*% S
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

P <- X0 %*% t(Y0)
min_p <- min(P)
max_p <- max(P)
tau <- 0.005

cat(sprintf("Edge Probability Range: [%.4f, %.4f]\n", min_p, max_p))
cat(sprintf("Selected tau: %.6f\n", tau))
cat(sprintf("Adaptive line search: ls_beta = 0.8 (iter 1-8) -> 0.4 (iter 9+)\n\n"))

# Storage for single replication result
results <- list(
  rep = rep_id,
  SSE = NA,
  fit_time = NA,
  coverage_true = NA,
  coverage_plugin = NA
)

rep_start <- Sys.time()
cat(sprintf("============== REPLICATION %d ==============\n", rep_id))

# Generate data
cat("Generating adjacency matrix A...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

cat("Generating covariate matrix B...\n")
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# Fit model with adaptive step size
cat("Fitting GRDPG model with ADAPTIVE step sizes...\n")
fit_start <- Sys.time()
fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                      maxit = maxit, tol = tol, tau = tau)
fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))

# Manual Oracle Procrustes alignment
X_est <- fit$X
M <- t(X_est) %*% X0
svd_res <- svd(M)
Q <- svd_res$u %*% t(svd_res$v)
X_aligned <- X_est %*% Q

Z_est <- fit$Z
Y_aligned <- X_aligned %*% S

# Refit Z with aligned X
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

SSE <- sum((X_aligned - X0)^2)

cat(sprintf("SSE: %.4f, Fit Time: %.1f seconds\n", SSE, fit_time))

# Compute coverage with BOTH methods on ALL nodes (parallelized)
cat(sprintf("Computing coverage on ALL nodes using %d cores...\n", n_cores))
chi2_crit <- qchisq(0.95, df = d)

# Function to compute coverage for a single node
compute_node_coverage <- function(i) {
  result <- list(true = NA, plugin = NA)

  # TRUE G_in
  G_in_true <- compute_G_in_true(i, X0, Y0, Z0, tau = tau)
  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
    result$true <- (mahal_dist <= chi2_crit)
  }

  # PLUG-IN G_in
  G_in_plugin <- compute_G_in_plugin(i, X_aligned, Y_aligned, Z_updated, tau = tau)
  eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_plugin %*% diff_vec)
    result$plugin <- (mahal_dist <= chi2_crit)
  }

  return(result)
}

# Parallel computation across all nodes
coverage_start <- Sys.time()
coverage_results <- mclapply(1:n, compute_node_coverage, mc.cores = n_cores)
coverage_time <- as.numeric(difftime(Sys.time(), coverage_start, units = "secs"))

cat(sprintf("Coverage computation completed in %.1f seconds\n", coverage_time))

# Extract results
in_ellipse_true <- sapply(coverage_results, function(x) x$true)
in_ellipse_plugin <- sapply(coverage_results, function(x) x$plugin)

coverage_true <- mean(in_ellipse_true, na.rm = TRUE)
coverage_plugin <- mean(in_ellipse_plugin, na.rm = TRUE)

# Store results
results$SSE <- SSE
results$fit_time <- fit_time
results$coverage_true <- coverage_true
results$coverage_plugin <- coverage_plugin

cat(sprintf("Coverage: TRUE=%.1f%%, PLUGIN=%.1f%%\n\n",
            100*coverage_true, 100*coverage_plugin))

rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

cat("============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETED\n", rep_id))
cat(sprintf("  Total time: %.2f minutes (Fit: %.1fs, Coverage: %.1fs)\n",
            rep_time, fit_time, coverage_time))
cat("============================================================================\n\n")

# Print results
cat("RESULTS:\n")
cat(sprintf("%-5s %-10s %-10s %-15s %-15s\n", "Rep", "SSE", "Time(s)", "Cov(TRUE)%", "Cov(PLUGIN)%"))
cat(paste(rep("-", 60), collapse=""), "\n")
cat(sprintf("%-5d %-10.4f %-10.1f %-15.1f %-15.1f\n",
            results$rep, results$SSE, results$fit_time,
            100*results$coverage_true, 100*results$coverage_plugin))

# Save results with replication number in filename
output_file <- sprintf("diagnostic_2d_tight_oracle_rep%d.rds", rep_id)
saveRDS(results, output_file)
cat(sprintf("\nResults saved to: %s\n", output_file))
cat("============================================================================\n")
