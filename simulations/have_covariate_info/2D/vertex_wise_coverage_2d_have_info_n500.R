#!/usr/bin/env Rscript
# Vertex-wise Coverage: cgrdpg vs ASE vs OSE with COVARIATE INFORMATION
# Scenario: Z0 = rnorm, B contains signal
# This tests cgrdpg performance when covariates are informative
# n=500, p_cov=250, d=2 (RDPG)

library(cgrdpg)

# --- Parameters ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript vertex_wise_coverage_2d_have_info_n500.R <rep_number>")
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) stop("Rep number must be between 1 and 100")

n         <- 500
p_cov     <- 250
d         <- 2
maxit     <- 30
tol       <- 0.01
tau       <- 0.005
eps_clip  <- 1e-10
chi2_crit <- qchisq(0.95, df = d)
S         <- diag(c(1, 1))  # RDPG (positive definite)

output_dir <- "results_2d_have_info_n500"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("============================================================================\n")
cat("  Vertex-wise Coverage: WITH COVARIATE INFORMATION (n=500)\n")
cat(sprintf("  Replication %d/100\n", rep_id))
cat("  Z0 = rnorm, B contains signal\n")
cat("  Methods: cgrdpg, ASE, OSE\n")
cat("============================================================================\n\n")

# --- Helpers ---

procrustes_align <- function(X_est, X_target) {
  svd_res <- svd(t(X_est) %*% X_target)
  Q <- svd_res$u %*% t(svd_res$v)
  list(X_aligned = X_est %*% Q, Q = Q)
}

compute_ose_step <- function(A, X_init, clipping_val) {
  n_nodes <- nrow(A)
  d_dim   <- ncol(X_init)
  X_new   <- matrix(0, n_nodes, d_dim)
  for (i in 1:n_nodes) {
    x_i <- X_init[i, ]
    idx_j <- setdiff(1:n_nodes, i)
    p_i <- pmax(pmin(as.vector(X_init[idx_j, ] %*% x_i), 1 - clipping_val), clipping_val)
    resid <- A[i, idx_j] - p_i
    w_score <- 1 / (p_i * (1 - p_i))
    grad <- colSums(X_init[idx_j, ] * (resid * w_score))
    G <- t(X_init[idx_j, ]) %*% (X_init[idx_j, ] * w_score)
    X_new[i, ] <- x_i + solve(G + diag(1e-9, d_dim), grad)
  }
  X_new
}

compute_G_in_cgrdpg <- function(i, X_mat, Y_mat, Z_mat, tau) {
  n_loc <- nrow(X_mat); p_loc <- nrow(Z_mat)
  s <- as.vector(X_mat[i, ] %*% t(Y_mat))
  w <- dpsi(s, tau = tau); w[i] <- 0
  G_net <- crossprod(Y_mat * sqrt(w))
  (G_net + crossprod(Z_mat)) / (n_loc + p_loc)
}

compute_prec_ase <- function(i, X_mat, S, clipping_val) {
  idx_j  <- setdiff(1:nrow(X_mat), i)
  p_vals <- pmax(pmin(as.vector(X_mat[idx_j, ] %*% (S %*% X_mat[i, ])), 1 - clipping_val), clipping_val)
  Delta  <- t(X_mat) %*% X_mat
  M_mat  <- t(X_mat[idx_j, ]) %*% (X_mat[idx_j, ] * p_vals * (1 - p_vals))
  S %*% Delta %*% solve(M_mat + diag(1e-9, d), Delta) %*% S
}

compute_prec_ose <- function(i, X_mat, clipping_val) {
  idx_j  <- setdiff(1:nrow(X_mat), i)
  p_vals <- pmax(pmin(as.vector(X_mat[idx_j, ] %*% X_mat[i, ]), 1 - clipping_val), clipping_val)
  w      <- 1 / (p_vals * (1 - p_vals))
  t(X_mat[idx_j, ]) %*% (X_mat[idx_j, ] * w)
}

check_coverage <- function(err, Prec, scale = 1.0) {
  ev <- eigen(Prec, only.values = TRUE)$values
  if (min(ev) < 1e-10) return(NA)
  (scale * as.numeric(t(err) %*% Prec %*% err)) <= chi2_crit
}

# ============================================================================
#  SINGLE REPLICATION
# ============================================================================
rep_start <- Sys.time()
set.seed(598 + rep_id)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (ncores <= 1) ncores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Using %d cores for cgrdpg parallel fitting\n\n", ncores))

  # 1. Latent positions (same as no_covariate sims)
  theta <- pi * (1:n) / (n - 1)
  X0 <- matrix(0, n, d)
  X0[, 1] <- 0.28 * sin(theta) + 0.42
  X0[, 2] <- 0.28 * cos(theta) + 0.42
  Y0 <- X0 %*% S

  # 2. WITH COVARIATE INFORMATION: Z0 = rnorm
  Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

  P  <- X0 %*% t(Y0)
  cat(sprintf("Edge probability range: [%.4f, %.4f]\n\n", min(P), max(P)))

  # 3. Data
  A <- (runif(n^2) < P) * 1.0
  A <- A * upper.tri(A, diag = FALSE) + t(A * upper.tri(A, diag = FALSE))

  # B contains SIGNAL (Z0 %*% t(X0)) + noise
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

  # 4. ASE
  cat("Computing ASE...\n")
  t0        <- Sys.time()
  A_aug     <- A; diag(A_aug) <- rowSums(A) / (n - 1)
  ase_fit   <- ase_grdpg(A_aug, d = d)
  X_ase_raw <- ase_fit$X
  ase_time  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  X_ase     <- procrustes_align(X_ase_raw, X0)$X_aligned
  cat(sprintf("ASE: time=%.1fs\n", ase_time))

  # 5. cgrdpg
  cat("Fitting cgrdpg...\n")
  t0  <- Sys.time()
  fit <- tryCatch(
    fit_grdpg_cov_parallel(A, B, d = d, p = 2, q = 0,
                           maxit = maxit, tol = tol, tau = tau, ncores = ncores),
    error = function(e) fit_grdpg_cov(A, B, d = d, p = 2, q = 0,
                                      maxit = maxit, tol = tol, tau = tau)
  )
  cgrdpg_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  X_cgrdpg    <- procrustes_align(fit$X, X0)$X_aligned
  Y_cgrdpg    <- X_cgrdpg %*% S
  Z_cgrdpg    <- fit$Z
  cat(sprintf("cgrdpg: converged=%s, iters=%d, time=%.1fs\n",
              fit$converged, fit$iters, cgrdpg_time))

  # 6. OSE
  cat("Computing OSE...\n")
  t0        <- Sys.time()
  X_ose_raw <- compute_ose_step(A, X_ase_raw, eps_clip)
  ose_step_time  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ose_time  <- ase_time + ose_step_time
  X_ose     <- procrustes_align(X_ose_raw, X0)$X_aligned
  cat(sprintf("OSE: time=%.1fs (ASE: %.1fs + step: %.1fs)\n", ose_time, ase_time, ose_step_time))

  sse <- c(cgrdpg = sum((X_cgrdpg - X0)^2),
           ase    = sum((X_ase    - X0)^2),
           ose    = sum((X_ose    - X0)^2))
  cat(sprintf("SSE  cgrdpg=%.4f  ASE=%.4f  OSE=%.4f\n\n", sse["cgrdpg"], sse["ase"], sse["ose"]))

  # 7. Vertex-wise coverage
  cat(sprintf("Computing vertex-wise coverage for all %d vertices...\n", n))
  results_mat <- matrix(NA_real_, nrow = n, ncol = 6,
    dimnames = list(NULL, c("cgrdpg_true", "cgrdpg_plugin",
                             "ase_true",    "ase_plugin",
                             "ose_true",    "ose_plugin")))

  t0 <- Sys.time()
  for (i in 1:n) {
    if (i %% 100 == 0) cat(sprintf("  Vertex %d/%d\n", i, n))

    # cgrdpg
    results_mat[i, "cgrdpg_true"]   <- check_coverage(
      X0[i,] - X_cgrdpg[i,],
      compute_G_in_cgrdpg(i, X0,       Y0,       Z0,       tau), n + p_cov)
    results_mat[i, "cgrdpg_plugin"] <- check_coverage(
      X0[i,] - X_cgrdpg[i,],
      compute_G_in_cgrdpg(i, X_cgrdpg, Y_cgrdpg, Z_cgrdpg, tau), n + p_cov)

    # ASE
    results_mat[i, "ase_true"]   <- check_coverage(
      X_ase[i,] - X0[i,], compute_prec_ase(i, X0,    S, eps_clip))
    results_mat[i, "ase_plugin"] <- check_coverage(
      X_ase[i,] - X0[i,], compute_prec_ase(i, X_ase, S, eps_clip))

    # OSE
    results_mat[i, "ose_true"]   <- check_coverage(
      X_ose[i,] - X0[i,], compute_prec_ose(i, X0,   eps_clip))
    results_mat[i, "ose_plugin"] <- check_coverage(
      X_ose[i,] - X0[i,], compute_prec_ose(i, X_ose, eps_clip))
  }

  cov_time    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  overall_cov <- colMeans(results_mat, na.rm = TRUE)
  rep_time    <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

  cat("\nOverall coverage (this rep):\n")
  for (nm in names(overall_cov))
    cat(sprintf("  %-20s %.1f%%  (NAs: %d)\n", nm, 100 * overall_cov[nm],
                sum(is.na(results_mat[, nm]))))
  cat(sprintf("\nTotal rep time: %.2f min\n", rep_time))

  # 8. Save
  out_file <- file.path(output_dir, sprintf("rep_%03d.rds", rep_id))
  saveRDS(list(
    rep_id      = rep_id,
    seed        = 598 + rep_id,
    n = n, p_cov = p_cov, d = d, tau = tau, S = S,
    results_mat = results_mat,
    overall_cov = overall_cov,
    n_na        = colSums(is.na(results_mat)),
    sse         = sse,
    timing      = list(cgrdpg_time = cgrdpg_time, ase_time = ase_time,
                       ose_time = ose_time, cov_time = cov_time,
                       rep_time_min = rep_time),
    converged  = fit$converged,
    iterations = fit$iters
  ), out_file)

cat(sprintf("\nResults saved to: %s\n", out_file))
cat("============================================================================\n")
