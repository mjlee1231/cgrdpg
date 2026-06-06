#!/usr/bin/env Rscript
# ASE vs OSE vs cgrdpg vertex-wise coverage: 3D GRDPG, n=1000, 100 reps
# Single node: all 100 replications run sequentially
# Methods: cgrdpg-TRUE, cgrdpg-PLUGIN, ASE-TRUE, ASE-PLUGIN, OSE-TRUE, OSE-PLUGIN
# Signature matrix: S = diag(1, 1, -1), p_cov=500

library(cgrdpg)

# --- Parameters ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript ase_ose_vertex_wise_3d_n1000.R <rep_number>")
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) stop("Rep number must be between 1 and 100")

n         <- 1000
p_cov     <- 500
d         <- 3
maxit     <- 30
tol       <- 0.01
tau       <- 0.001
eps_clip  <- 1e-10
chi2_crit <- qchisq(0.95, df = d)
S         <- diag(c(1, 1, -1))

output_dir <- "results_3d_ase_ose_cgrdpg_n1000"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("============================================================================\n")
cat("  ASE vs OSE vs cgrdpg Vertex-wise Coverage: 3D GRDPG, n=1000\n")
cat(sprintf("  Replication %d/100\n", rep_id))
cat("  Methods: cgrdpg-TRUE, cgrdpg-PLUGIN, ASE-TRUE, ASE-PLUGIN, OSE-TRUE, OSE-PLUGIN\n")
cat("  S = diag(1, 1, -1), p_cov=500\n")
cat("============================================================================\n\n")

# --- Helpers ---

procrustes_align <- function(X_est, X_target) {
  svd_res <- svd(t(X_est) %*% X_target)
  Q <- svd_res$u %*% t(svd_res$v)
  list(X_aligned = X_est %*% Q, Q = Q)
}

compute_ose_step <- function(A, X_unsigned, X_signed, clipping_val) {
  # OSE for GRDPG
  # Edge prob: p_ji = x̃_j^T * x̂_i (signed × unsigned)
  # Gradient & Fisher info use signed ASE
  # Update: x_new = x_unsigned + G^{-1} * grad (PLUS sign)
  n_nodes <- nrow(A)
  d_dim   <- ncol(X_unsigned)
  X_new   <- matrix(0, n_nodes, d_dim)
  for (i in 1:n_nodes) {
    x_i_unsigned <- X_unsigned[i, ]
    idx_j   <- setdiff(1:n_nodes, i)
    # Edge probability: p_ji = x̃_j^T * x̂_i
    p_i     <- pmax(pmin(as.vector(X_signed[idx_j, ] %*% x_i_unsigned), 1 - clipping_val), clipping_val)
    resid   <- A[i, idx_j] - p_i
    w_score <- 1 / (p_i * (1 - p_i))
    # Gradient using signed ASE
    grad    <- colSums(X_signed[idx_j, ] * (resid * w_score))
    # Fisher information using signed ASE
    G       <- t(X_signed[idx_j, ]) %*% (X_signed[idx_j, ] * w_score)
    # Newton-Raphson update from unsigned ASE
    X_new[i, ] <- x_i_unsigned + solve(G + diag(1e-9, d_dim), grad)
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
  # p_ij = x_i^T S x_j
  p_vals <- pmax(pmin(as.vector(X_mat[idx_j, ] %*% (S %*% X_mat[i, ])), 1 - clipping_val), clipping_val)
  # Delta = X^T X (second moment matrix)
  Delta  <- t(X_mat) %*% X_mat
  # M = X_{-i}^T diag(p*(1-p)) X_{-i}  (sandwich middle)
  M_mat  <- t(X_mat[idx_j, ]) %*% (X_mat[idx_j, ] * p_vals * (1 - p_vals))
  # ASE precision: S %*% Delta %*% M^{-1} %*% Delta %*% S
  S %*% Delta %*% solve(M_mat + diag(1e-9, d), Delta) %*% S
}

compute_prec_ose <- function(i, X_unsigned, X_signed, clipping_val) {
  # OSE precision for GRDPG
  # Edge prob: p_ji = x̃_j^T * x̂_i (signed × unsigned)
  # Precision: X̃^T * diag(w) * X̃ (using signed ASE)
  idx_j  <- setdiff(1:nrow(X_unsigned), i)
  # Edge probability: p_ji = x̃_j^T * x̂_i
  p_vals <- pmax(pmin(as.vector(X_signed[idx_j, ] %*% X_unsigned[i, ]), 1 - clipping_val), clipping_val)
  w      <- 1 / (p_vals * (1 - p_vals))
  # Precision using signed ASE
  t(X_signed[idx_j, ]) %*% (X_signed[idx_j, ] * w)
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

  # 1. Latent positions
  t <- (0:(n - 1)) / n
  X0 <- cbind(0.15 * sin(2*pi*t) + 0.6,
              0.15 * cos(2*pi*t) + 0.6,
              0.15 * cos(4*pi*t))
  Y0 <- X0 %*% S
  Z0 <- matrix(rnorm(p_cov * d), p_cov, d)
  P  <- X0 %*% t(Y0)
  cat(sprintf("Edge probability range: [%.4f, %.4f]\n\n", min(P), max(P)))

  # 2. Data
  A <- (runif(n^2) < P) * 1.0
  A <- A * upper.tri(A, diag = FALSE) + t(A * upper.tri(A, diag = FALSE))
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

  # 3. ASE (computed first to get estimated signature matrix)
  cat("Computing ASE...\n")
  t0        <- Sys.time()
  A_aug     <- A; diag(A_aug) <- rowSums(A) / (n - 1)
  ase_fit   <- ase_grdpg(A_aug, d = d)
  X_ase_unsigned <- ase_fit$X         # U |Λ|^{1/2}
  X_ase_signed   <- ase_fit$X_signed  # U |Λ|^{1/2} sign(Λ)
  S_estimated    <- ase_fit$sign_diag # Signature from actual eigenvalues
  ase_time  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  X_ase     <- procrustes_align(X_ase_unsigned, X0)$X_aligned
  cat(sprintf("ASE: time=%.1fs\n", ase_time))

  # 4. cgrdpg
  cat("Fitting cgrdpg...\n")
  t0  <- Sys.time()
  fit <- tryCatch(
    fit_grdpg_cov_parallel(A, B, d = d, p = 2, q = 1,
                           maxit = maxit, tol = tol, tau = tau, ncores = ncores),
    error = function(e) fit_grdpg_cov(A, B, d = d, p = 2, q = 1,
                                      maxit = maxit, tol = tol, tau = tau)
  )
  cgrdpg_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  X_cgrdpg    <- procrustes_align(fit$X, X0)$X_aligned
  Y_cgrdpg    <- X_cgrdpg %*% S_estimated  # Use estimated signature
  Z_cgrdpg    <- B %*% X_cgrdpg %*% solve(t(X_cgrdpg) %*% X_cgrdpg)
  cat(sprintf("cgrdpg: converged=%s, iters=%d, time=%.1fs\n",
              fit$converged, fit$iters, cgrdpg_time))

  # 5. OSE
  cat("Computing OSE...\n")
  t0        <- Sys.time()
  X_ose_raw <- compute_ose_step(A, X_ase_unsigned, X_ase_signed, eps_clip)
  ose_time  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  X_ose     <- procrustes_align(X_ose_raw, X0)$X_aligned
  cat(sprintf("OSE: time=%.1fs\n", ose_time))

  sse <- c(cgrdpg = sum((X_cgrdpg - X0)^2),
           ase    = sum((X_ase    - X0)^2),
           ose    = sum((X_ose    - X0)^2))
  cat(sprintf("SSE  cgrdpg=%.4f  ASE=%.4f  OSE=%.4f\n\n", sse["cgrdpg"], sse["ase"], sse["ose"]))

  # 6. Vertex-wise coverage
  cat(sprintf("Computing vertex-wise coverage for all %d vertices...\n", n))
  results_mat <- matrix(NA_real_, nrow = n, ncol = 6,
    dimnames = list(NULL, c("cgrdpg_true", "cgrdpg_plugin",
                             "ase_true",    "ase_plugin",
                             "ose_true",    "ose_plugin")))

  t0 <- Sys.time()
  for (i in 1:n) {
    if (i %% 200 == 0) cat(sprintf("  Vertex %d/%d\n", i, n))

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

    # OSE (use unsigned and signed versions)
    results_mat[i, "ose_true"]   <- check_coverage(
      X_ose[i,] - X0[i,], compute_prec_ose(i, X0, Y0, eps_clip))
    results_mat[i, "ose_plugin"] <- check_coverage(
      X_ose[i,] - X0[i,], compute_prec_ose(i, X_ose, X_ose %*% S_estimated, eps_clip))
  }

  cov_time    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  overall_cov <- colMeans(results_mat, na.rm = TRUE)
  rep_time    <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

  cat("\nOverall coverage (this rep):\n")
  for (nm in names(overall_cov))
    cat(sprintf("  %-20s %.1f%%  (NAs: %d)\n", nm, 100 * overall_cov[nm],
                sum(is.na(results_mat[, nm]))))
  cat(sprintf("\nTotal rep time: %.2f min\n", rep_time))

  # 7. Save
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
