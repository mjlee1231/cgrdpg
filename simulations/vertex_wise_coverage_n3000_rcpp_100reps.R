#!/usr/bin/env Rscript
# Vertex-wise Coverage with PARALLEL + RCPP optimization: n=3000, 100 reps
# Uses fit_grdpg_cov_parallel for 4x parallel speedup
# Uses Rcpp functions (dpsi_cpp, compute_G_in_vectorized_cpp) for 2x additional speedup
# Total expected speedup: ~6-8x vs serial
# Expected time per rep: ~40-50 minutes (vs ~180 min with serial)

library(cgrdpg)

# Get replication number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript vertex_wise_coverage_n3000_rcpp_100reps.R <rep_number>")
}
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) {
  stop("Rep number must be between 1 and 100")
}

# Check if Rcpp functions are available
rcpp_available <- exists("dpsi_cpp") && exists("compute_G_in_vectorized_cpp")
if (!rcpp_available) {
  cat("WARNING: Rcpp functions not found! Falling back to R versions.\n")
  cat("Performance will be degraded. Please reinstall package with Rcpp support.\n\n")
}

cat("============================================================================\n")
cat("  VERTEX-WISE COVERAGE with PARALLEL + RCPP: n=3000\n")
cat(sprintf("  Running REPLICATION %d/100\n", rep_id))
cat("  Optimization: Parallel fitting + Rcpp C++ functions\n")
cat(sprintf("  Rcpp enabled: %s\n", ifelse(rcpp_available, "YES", "NO (DEGRADED)")))
cat("============================================================================\n\n")

# Parameters
n <- 3000
p_cov <- 1500
d <- 2
maxit <- 30
tol <- 0.01
tau <- 0.005

# Detect available cores (will be set by SLURM)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (ncores <= 1) {
  ncores <- parallel::detectCores() - 1
}
ncores <- max(1, ncores)  # At least 1 core

cat(sprintf("Available cores: %d\n", ncores))

set.seed(598 + rep_id)

# Generate data
cat("Generating 2D SPREAD latent positions (n=3000)...\n")
i_vals <- 1:n
theta <- pi * i_vals / (n - 1)
X0 <- matrix(0, n, 2)
X0[, 1] <- 0.28 * sin(theta) + 0.42
X0[, 2] <- 0.28 * cos(theta) + 0.42

S <- diag(c(1, 1))
Y0 <- X0 %*% S
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

P <- X0 %*% t(Y0)
cat(sprintf("Edge Probability Range: [%.4f, %.4f]\n", min(P), max(P)))
cat(sprintf("Parameters: n=%d, p_cov=%d, d=%d, tau=%.3f\n", n, p_cov, d, tau))
cat(sprintf("Cores for parallel fitting: %d\n\n", ncores))

rep_start <- Sys.time()

# Generate adjacency and covariate matrices
cat("Generating adjacency matrix A...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

cat("Generating covariate matrix B...\n")
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# Fit model using PARALLEL version
cat(sprintf("Fitting GRDPG model (PARALLEL with %d cores)...\n", ncores))
fit_start <- Sys.time()

fit <- fit_grdpg_cov_parallel(A, B, d = d, p = 2, q = 0,
                               maxit = maxit, tol = tol, tau = tau,
                               ncores = ncores)

fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
cat(sprintf("Model fit completed: Converged=%s, Iterations=%d, Time=%.1f sec\n",
            fit$converged, fit$iters, fit_time))

# Oracle Procrustes alignment
cat("Performing Procrustes alignment...\n")
X_est <- fit$X
M <- t(X_est) %*% X0
svd_res <- svd(M)
Q <- svd_res$u %*% t(svd_res$v)
X_aligned <- X_est %*% Q

Z_est <- fit$Z
Y_aligned <- X_aligned %*% S
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

SSE <- sum((X_aligned - X0)^2)
cat(sprintf("SSE: %.4f\n\n", SSE))

# Compute coverage for ALL 3000 vertices using RCPP if available
cat("Computing elliptical coverage for ALL 3000 vertices...\n")
cat(sprintf("Using %s G_in computation\n",
            ifelse(rcpp_available, "Rcpp (C++)", "R (slower)")))
chi2_crit <- qchisq(0.95, df = d)

in_ellipse_true <- logical(n)
in_ellipse_plugin <- logical(n)

coverage_start <- Sys.time()

for (i in 1:n) {
  if (i %% 500 == 0) cat(sprintf("  Vertex %d/%d (%.1f%%)\n", i, n, 100*i/n))

  # TRUE G_in - use Rcpp if available
  if (rcpp_available) {
    G_in_true <- compute_G_in_vectorized_cpp(i - 1, X0, Y0, Z0, tau = tau)
  } else {
    # Fallback to R version
    s <- as.vector(X0[i, ] %*% t(Y0))
    w <- dpsi(s, tau = tau)
    w[i] <- 0
    Y_weighted <- Y0 * sqrt(w)
    G_net <- crossprod(Y_weighted)
    G_cov <- crossprod(Z0)
    G_in_true <- (G_net + G_cov) / (n + p_cov)
  }

  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) >= 1e-10) {
    diff_vec <- X0[i, ] - X_aligned[i, ]
    mahal_dist <- (n + p_cov) * as.numeric(t(diff_vec) %*% G_in_true %*% diff_vec)
    in_ellipse_true[i] <- (mahal_dist <= chi2_crit)
  } else {
    in_ellipse_true[i] <- NA
  }

  # PLUG-IN G_in - use Rcpp if available
  if (rcpp_available) {
    G_in_plugin <- compute_G_in_vectorized_cpp(i - 1, X_aligned, Y_aligned, Z_updated, tau = tau)
  } else {
    # Fallback to R version
    s <- as.vector(X_aligned[i, ] %*% t(Y_aligned))
    w <- dpsi(s, tau = tau)
    w[i] <- 0
    Y_weighted <- Y_aligned * sqrt(w)
    G_net <- crossprod(Y_weighted)
    G_cov <- crossprod(Z_updated)
    G_in_plugin <- (G_net + G_cov) / (n + p_cov)
  }

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

# Compute overall coverage
coverage_true <- mean(in_ellipse_true, na.rm = TRUE)
coverage_plugin <- mean(in_ellipse_plugin, na.rm = TRUE)
n_na_true <- sum(is.na(in_ellipse_true))
n_na_plugin <- sum(is.na(in_ellipse_plugin))

cat(sprintf("Coverage (this rep): TRUE=%.1f%% (NAs=%d), PLUGIN=%.1f%% (NAs=%d)\n",
            100*coverage_true, n_na_true, 100*coverage_plugin, n_na_plugin))

rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

# Save results
output_dir <- "vertex_wise_results_n3000_rcpp"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

results <- list(
  rep = rep_id,
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
  ncores = ncores,
  rcpp_enabled = rcpp_available,
  in_ellipse_true = in_ellipse_true,
  in_ellipse_plugin = in_ellipse_plugin,
  SSE = SSE,
  fit_time = fit_time,
  coverage_time = coverage_time,
  coverage_true = coverage_true,
  coverage_plugin = coverage_plugin,
  n_na_true = n_na_true,
  n_na_plugin = n_na_plugin,
  parallel = TRUE
)

output_file <- file.path(output_dir, sprintf("vertex_wise_n3000_rep%03d.rds", rep_id))
saveRDS(results, output_file)

cat("\n============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETED\n", rep_id))
cat(sprintf("  Total time: %.2f minutes (Fit: %.1fs, Coverage: %.1fs)\n",
            rep_time, fit_time, coverage_time))
cat(sprintf("  Optimization: Parallel (%d cores) + Rcpp C++\n", ncores))
cat(sprintf("  Expected speedup: ~6-8x vs serial\n"))
cat(sprintf("  Results saved to: %s\n", output_file))
cat("============================================================================\n")
