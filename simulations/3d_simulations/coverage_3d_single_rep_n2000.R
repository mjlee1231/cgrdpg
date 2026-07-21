#!/usr/bin/env Rscript
# coverage_3d_single_rep_n2000.R
# Single replication for 3D GRDPG coverage using R package (n=2000, p_cov=1000)
# Uses OPTIMIZED latent position design: X0 = [0.42*t+0.46, 0.27*sin+0.46, 0.20*cos]

library(cgrdpg)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript coverage_3d_single_rep_n2000.R <rep_number>")
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) stop("Rep number must be between 1 and 100")

# Parameters
n         <- 2000
p_cov     <- 1000
d         <- 3
maxit     <- 30
tol       <- 0.01
tau       <- 0.001
eps_clip  <- 1e-10
chi2_crit <- qchisq(0.95, df = d)
S         <- diag(c(1, 1, -1))

output_dir <- "results_r_3d_coverage_n2000"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("============================================================================\n")
cat("  cgrdpg Coverage Simulation (R): n=2000, p_cov=1000\n")
cat(sprintf("  Replication %d/100\n", rep_id))
cat("  Using OPTIMIZED latent position design\n")
cat("  S = diag(1, 1, -1)\n")
cat("============================================================================\n\n")

# --- Helper Functions ---

procrustes_align <- function(X_est, X_target) {
  svd_res <- svd(t(X_est) %*% X_target)
  Q <- svd_res$u %*% t(svd_res$v)
  list(X_aligned = X_est %*% Q, Q = Q)
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

# 1. OPTIMIZED Latent positions (λ3=20, condition#=33.3)
t <- (1:n) / n
X0 <- cbind(
  0.42 * t + 0.46,
  0.27 * sin(2*pi*t) + 0.46,
  0.20 * cos(4*pi*t)
)
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
X_ase_unsigned <- ase_fit$X
X_ase_signed   <- ase_fit$X_signed
S_estimated    <- ase_fit$sign_diag
ase_time  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
X_ase     <- procrustes_align(X_ase_unsigned, X0)$X_aligned
cat(sprintf("ASE: time=%.1fs, S_est=diag([%+d,%+d,%+d])\n",
            ase_time, S_estimated[1,1], S_estimated[2,2], S_estimated[3,3]))

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

# SSE
sse_cgrdpg <- sum((X_cgrdpg - X0)^2)
sse_ase    <- sum((X_ase - X0)^2)
cat(sprintf("SSE: cgrdpg=%.4f  ASE=%.4f\n\n", sse_cgrdpg, sse_ase))

# 5. Vertex-wise coverage
cat(sprintf("Computing vertex-wise coverage for all %d vertices...\n", n))
results_mat <- matrix(NA_real_, nrow = n, ncol = 4,
  dimnames = list(NULL, c("cgrdpg_true", "cgrdpg_plugin",
                           "ase_true",    "ase_plugin")))

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
    X_ase[i,] - X0[i,], compute_prec_ase(i, X_ase, S_estimated, eps_clip))
}
cov_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Overall coverage
overall_cov <- colMeans(results_mat, na.rm = TRUE)
n_na <- colSums(is.na(results_mat))

rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "mins"))

cat("\nOverall coverage (this rep):\n")
cat(sprintf("  cgrdpg-TRUE:    %.1f%%  (NAs: %d)\n", 100 * overall_cov["cgrdpg_true"], n_na["cgrdpg_true"]))
cat(sprintf("  cgrdpg-PLUGIN:  %.1f%%  (NAs: %d)\n", 100 * overall_cov["cgrdpg_plugin"], n_na["cgrdpg_plugin"]))
cat(sprintf("  ASE-TRUE:       %.1f%%  (NAs: %d)\n", 100 * overall_cov["ase_true"], n_na["ase_true"]))
cat(sprintf("  ASE-PLUGIN:     %.1f%%  (NAs: %d)\n", 100 * overall_cov["ase_plugin"], n_na["ase_plugin"]))
cat(sprintf("\nTotal rep time: %.2f min\n", rep_time))

# Save results
results <- list(
  rep_id = rep_id,
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
  S = S,
  S_estimated = S_estimated,
  results_mat = results_mat,
  overall_cov = overall_cov,
  n_na = n_na,
  sse_cgrdpg = sse_cgrdpg,
  sse_ase = sse_ase,
  cgrdpg_time = cgrdpg_time,
  ase_time = ase_time,
  cov_time = cov_time,
  rep_time = rep_time,
  X0 = X0,
  X_cgrdpg = X_cgrdpg,
  Y_cgrdpg = Y_cgrdpg,
  Z_cgrdpg = Z_cgrdpg,
  X_ase = X_ase,
  converged = fit$converged,
  iters = fit$iters
)

saveRDS(results, file.path(output_dir, sprintf("rep_%03d.rds", rep_id)))

cat(sprintf("\nResults saved to: %s\n",
            file.path(output_dir, sprintf("rep_%03d.rds", rep_id))))
cat("============================================================================\n")
