#!/usr/bin/env Rscript
# Vertex-wise Coverage Diagnostic: 2D SPREAD with Oracle Procrustes
# Computes coverage for ALL 500 vertices - saves per-vertex indicators
# To be aggregated across 100 replications for vertex-wise coverage rates
# ARRAY JOB VERSION: Each job runs 1 replication
#
# Configuration: Optimized SPREAD (r=0.28, c=0.42)
# Edge probability range: [0.195, 0.764]

library(cgrdpg)

# Get replication number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript vertex_wise_coverage_spread_100reps.R <rep_number>")
}
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) {
  stop("Rep number must be between 1 and 100")
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
cat("  VERTEX-WISE COVERAGE: 2D SPREAD - Oracle Procrustes\n")
cat(sprintf("  Running REPLICATION %d/100\n", rep_id))
cat("  Computing coverage for ALL 500 vertices\n")
cat("============================================================================\n\n")

# Fixed parameters (as specified)
n <- 500
p_cov <- 250
d <- 2
maxit <- 30
tol <- 0.01
tau <- 0.005
# Note: ls_beta = 0.4 was requested but is not a parameter in fit_grdpg_cov
# The function uses its default line search behavior

# Set seed for this replication (base seed 598)
set.seed(598 + rep_id)

# Generate 2D latent positions (SPREAD configuration - OPTIMIZED)
# Target edge probability range: [0.195, 0.764]
cat("Generating 2D SPREAD latent positions (optimized)...\n")
i_vals <- 1:n
theta <- pi * i_vals / (n - 1)
X0 <- matrix(0, n, 2)
X0[, 1] <- 0.28 * sin(theta) + 0.42
X0[, 2] <- 0.28 * cos(theta) + 0.42

S <- diag(c(1, 1))
Y0 <- X0 %*% S
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

P <- X0 %*% t(Y0)
min_p <- min(P)
max_p <- max(P)

cat(sprintf("Edge Probability Range: [%.4f, %.4f]\n", min_p, max_p))
cat(sprintf("Parameters: n=%d, p_cov=%d, d=%d, tau=%.3f\n\n",
            n, p_cov, d, tau))

rep_start <- Sys.time()

# Generate data
cat("Generating adjacency matrix A...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

cat("Generating covariate matrix B...\n")
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# Fit model
cat("Fitting GRDPG model...\n")
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

# Compute coverage for ALL 500 vertices using ELLIPTICAL confidence intervals
cat("Computing elliptical coverage for ALL 500 vertices...\n")
chi2_crit <- qchisq(0.95, df = d)

# Initialize storage for ALL vertices
in_ellipse_true <- logical(n)
in_ellipse_plugin <- logical(n)

coverage_start <- Sys.time()

for (i in 1:n) {
  if (i %% 100 == 0) cat(sprintf("  Vertex %d/%d\n", i, n))

  # TRUE G_in (using true latent positions)
  G_in_true <- compute_G_in_true(i, X0, Y0, Z0, tau = tau)
  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
    in_ellipse_true[i] <- (mahal_dist <= chi2_crit)
  } else {
    in_ellipse_true[i] <- NA
  }

  # PLUG-IN G_in (using estimated latent positions)
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

coverage_time <- as.numeric(difftime(Sys.time(), coverage_start, units = "secs"))
cat(sprintf("Coverage computation completed in %.1f seconds\n", coverage_time))

# Compute overall coverage for this rep
coverage_true <- mean(in_ellipse_true, na.rm = TRUE)
coverage_plugin <- mean(in_ellipse_plugin, na.rm = TRUE)
n_na_true <- sum(is.na(in_ellipse_true))
n_na_plugin <- sum(is.na(in_ellipse_plugin))

cat(sprintf("Coverage (this rep): TRUE=%.1f%% (NAs=%d), PLUGIN=%.1f%% (NAs=%d)\n",
            100*coverage_true, n_na_true, 100*coverage_plugin, n_na_plugin))

rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

# Create output directory if it doesn't exist
# Using "_optimized" suffix to distinguish from original latent positions
output_dir <- "vertex_wise_results_n500_optimized"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save results with per-vertex indicators for ALL 500 vertices
results <- list(
  rep = rep_id,
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
  in_ellipse_true = in_ellipse_true,      # Length 500: TRUE/FALSE/NA for each vertex
  in_ellipse_plugin = in_ellipse_plugin,  # Length 500: TRUE/FALSE/NA for each vertex
  SSE = SSE,
  fit_time = fit_time,
  coverage_time = coverage_time,
  coverage_true = coverage_true,
  coverage_plugin = coverage_plugin,
  n_na_true = n_na_true,
  n_na_plugin = n_na_plugin
)

output_file <- file.path(output_dir, sprintf("vertex_wise_spread_rep%03d.rds", rep_id))
saveRDS(results, output_file)

cat("\n============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETED\n", rep_id))
cat(sprintf("  Total time: %.2f minutes (Fit: %.1fs, Coverage: %.1fs)\n",
            rep_time, fit_time, coverage_time))
cat(sprintf("  Results saved to: %s\n", output_file))
cat("============================================================================\n")
