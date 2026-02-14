#!/usr/bin/env Rscript
# ASE vs OSE vertex-wise coverage simulation: n=3000
# This script performs a single replication of a 2D GRDPG simulation.
# It computes vertex-wise coverage for ASE and OSE using manual Procrustes alignment.
# No covariates - pure network model

# Load necessary library
library(cgrdpg)

# --- 1. Simulation Parameters ---
args <- commandArgs(trailingOnly = TRUE)
rep_id <- if(length(args) > 0) as.numeric(args[1]) else 1
set.seed(598 + rep_id)  # Match seed with covariate simulations

n <- 3000
d <- 2
tau <- 0.005    # Convergence threshold for fit_grdpg_cov
eps_clip <- 1e-10 # Numerical clipping for probability weights
q_alpha <- qchisq(0.95, df = d) # 95% Confidence threshold for df=2

cat("============================================================================\n")
cat("  ASE vs OSE Vertex-wise Coverage Simulation\n")
cat(sprintf("  n=%d (no covariates)\n", n))
cat(sprintf("  Replication %d/100\n", rep_id))
cat(sprintf("  Seed: %d\n", 598 + rep_id))
cat("============================================================================\n\n")

# --- 2. Latent Position Generation (X0) ---
# Using optimized configuration: r=0.28, c=0.42 (same as covariate sims)
theta <- pi * (1:n) / (n-1)
X0 <- matrix(0, nrow = n, ncol = d)
X0[, 1] <- 0.28 * sin(theta) + 0.42
X0[, 2] <- 0.28 * cos(theta) + 0.42

# Generate Probability Matrix P (True model)
P <- X0 %*% t(X0)

cat(sprintf("Edge probability range: [%.4f, %.4f]\n", min(P), max(P)))
cat(sprintf("Parameters: n=%d, d=%d\n\n", n, d))

# Generate Observed Adjacency Matrix A (Bernoulli)
A <- matrix(runif(n^2), n, n) < P
A <- A * 1.0
diag(A) <- 0
A[lower.tri(A)] <- t(A)[lower.tri(A)]

# --- 3. OSE 1-Step Estimation Function ---
compute_ose_step <- function(A, X_init, clipping_val) {
  n_nodes <- nrow(A)
  d_dim <- ncol(X_init)
  X_new <- matrix(0, n_nodes, d_dim)

  for (i in 1:n_nodes) {
    x_i <- X_init[i, ]
    indices_j <- setdiff(1:n_nodes, i)

    # Calculate probabilities for the current node
    p_i <- as.vector(X_init[indices_j, ] %*% x_i)
    p_i <- pmax(pmin(p_i, 1 - clipping_val), clipping_val)

    # Score function (Gradient)
    resid <- A[i, indices_j] - p_i
    weight_score <- 1 / (p_i * (1 - p_i))
    grad <- colSums(X_init[indices_j, ] * (resid * weight_score))

    # Expected Fisher Information (Hessian G)
    G <- t(X_init[indices_j, ]) %*% (X_init[indices_j, ] * weight_score)

    # Newton-Raphson 1-step update
    step <- solve(G + diag(1e-9, d_dim), grad)
    X_new[i, ] <- x_i + step
  }
  return(X_new)
}

# --- 4. Adjacency Spectral Embedding (ASE) Estimation ---
cat("Computing ASE...\n")
ase_start <- Sys.time()

A_aug <- A
diag(A_aug) <- rowSums(A) / (n - 1)
ase_res <- ase_grdpg(A_aug, d = d)
X_ase_raw <- ase_res$X

ase_time <- as.numeric(difftime(Sys.time(), ase_start, units = "secs"))
cat(sprintf("ASE completed in %.1f seconds\n", ase_time))

# --- 5. OSE Estimation (1-Step Update) ---
cat("Computing OSE (1-step from ASE)...\n")
ose_start <- Sys.time()

X_ose_raw <- compute_ose_step(A, X_ase_raw, clipping_val = eps_clip)

ose_time <- as.numeric(difftime(Sys.time(), ose_start, units = "secs"))
cat(sprintf("OSE completed in %.1f seconds\n", ose_time))

# --- 6. Manual Oracle Procrustes Alignment ---
procrustes_align <- function(X_est, X_target) {
  M <- t(X_est) %*% X_target
  svd_res <- svd(M)
  Q <- svd_res$u %*% t(svd_res$v)
  return(X_est %*% Q)
}

X_ase_aligned <- procrustes_align(X_ase_raw, X0)
X_ose_aligned <- procrustes_align(X_ose_raw, X0)

# SSE Calculation
sse_ase <- sum((X_ase_aligned - X0)^2)
sse_ose <- sum((X_ose_aligned - X0)^2)

cat(sprintf("\nSSE: ASE=%.4f, OSE=%.4f\n\n", sse_ase, sse_ose))

# --- 7. Vertex-wise Coverage Calculation ---
cat(sprintf("Computing vertex-wise coverage for all %d vertices...\n", n))

get_precisions <- function(i, X_target_for_P) {
  idx_j <- setdiff(1:n, i)
  p_vals <- as.vector(X_target_for_P %*% X_target_for_P[i, ])
  p_vals <- pmax(pmin(p_vals, 1 - eps_clip), eps_clip)

  # OSE Precision (Fisher Information G)
  w_ose <- 1 / (p_vals[idx_j] * (1 - p_vals[idx_j]))
  prec_ose <- t(X_target_for_P[idx_j, ]) %*% (X_target_for_P[idx_j, ] * w_ose)

  # ASE Precision (Sandwich form: Delta %*% M^-1 %*% Delta)
  Delta <- t(X_target_for_P) %*% X_target_for_P
  w_ase <- p_vals[idx_j] * (1 - p_vals[idx_j])
  M_mat <- t(X_target_for_P[idx_j, ]) %*% (X_target_for_P[idx_j, ] * w_ase)
  prec_ase <- Delta %*% solve(M_mat + diag(1e-9, d)) %*% Delta

  return(list(prec_ose = prec_ose, prec_ase = prec_ase))
}

coverage_df <- data.frame(
  node = 1:n,
  ose_true = FALSE, ose_plugin = FALSE,
  ase_true = FALSE, ase_plugin = FALSE
)

coverage_start <- Sys.time()

for (i in 1:n) {
  if (i %% 500 == 0) cat(sprintf("  Vertex %d/%d (%.1f%%)\n", i, n, 100*i/n))

  # Calculate Precision Matrices
  p_true <- get_precisions(i, X0)
  p_ose_plug <- get_precisions(i, X_ose_aligned)
  p_ase_plug <- get_precisions(i, X_ase_aligned)

  # Calculate Error Vectors
  err_ose <- X_ose_aligned[i, ] - X0[i, ]
  err_ase <- X_ase_aligned[i, ] - X0[i, ]

  # Coverage Tests via Mahalanobis Distance (NO scaling by n)
  coverage_df$ose_true[i]   <- (t(err_ose) %*% p_true$prec_ose %*% err_ose) <= q_alpha
  coverage_df$ase_true[i]   <- (t(err_ase) %*% p_true$prec_ase %*% err_ase) <= q_alpha
  coverage_df$ose_plugin[i] <- (t(err_ose) %*% p_ose_plug$prec_ose %*% err_ose) <= q_alpha
  coverage_df$ase_plugin[i] <- (t(err_ase) %*% p_ase_plug$prec_ase %*% err_ase) <= q_alpha
}

coverage_time <- as.numeric(difftime(Sys.time(), coverage_start, units = "secs"))
cat(sprintf("Coverage computation completed in %.1f seconds\n", coverage_time))

# Compute overall coverage
cov_ase_true <- mean(coverage_df$ase_true)
cov_ase_plugin <- mean(coverage_df$ase_plugin)
cov_ose_true <- mean(coverage_df$ose_true)
cov_ose_plugin <- mean(coverage_df$ose_plugin)

cat(sprintf("\nOverall coverage (this rep):\n"))
cat(sprintf("  ASE-TRUE:   %.1f%%\n", 100*cov_ase_true))
cat(sprintf("  ASE-PLUGIN: %.1f%%\n", 100*cov_ase_plugin))
cat(sprintf("  OSE-TRUE:   %.1f%%\n", 100*cov_ose_true))
cat(sprintf("  OSE-PLUGIN: %.1f%%\n", 100*cov_ose_plugin))

total_time <- as.numeric(difftime(Sys.time(), ase_start, units = "secs"))

# --- 8. Export Results ---
# Create output directory if it doesn't exist
output_dir <- "ase_ose_results_n3000"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

result_bundle <- list(
  rep_id = rep_id,
  seed = 598 + rep_id,
  n = n,
  d = d,
  sse = list(ase = sse_ase, ose = sse_ose),
  coverage_df = coverage_df,
  overall_coverage = list(
    ase_true = cov_ase_true,
    ase_plugin = cov_ase_plugin,
    ose_true = cov_ose_true,
    ose_plugin = cov_ose_plugin
  ),
  timing = list(
    ase_time = ase_time,
    ose_time = ose_time,
    coverage_time = coverage_time,
    total_time = total_time
  )
)

output_file <- file.path(output_dir, sprintf("ase_ose_n3000_rep%03d.rds", rep_id))
saveRDS(result_bundle, output_file)

cat(sprintf("\nTotal time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))
cat(sprintf("Results saved to: %s\n", output_file))
cat("============================================================================\n")
