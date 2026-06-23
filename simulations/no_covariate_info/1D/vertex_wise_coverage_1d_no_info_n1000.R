#!/usr/bin/env Rscript
# 1D No Covariate Info: cgrdpg vs ASE vs OSE (SINGLE REPLICATION)
# Scenario: Z0 = 0, B = pure noise (no signal)
# Latent position: X_i = 0.75 * sin(π * i/(n-1)) + 0.1
# n=1000, p_cov=500

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript vertex_wise_coverage_1d_no_info_n1000.R <rep_id>")
}

rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) {
  stop("rep_id must be between 1 and 100")
}

library(cgrdpg)

cat("============================================================================\n")
cat(sprintf("  1D NO COVARIATE INFO - REPLICATION %d/100\n", rep_id))
cat("  Z0 = 0, B = pure noise\n")
cat("  Comparing: cgrdpg vs ASE vs OSE\n")
cat("============================================================================\n\n")

# Fixed parameters
n <- 1000
p_cov <- 500
d <- 1
maxit <- 30
tol <- 0.01
base_seed <- 598
eps_clip <- 1e-10
tau <- 0.005

# Set seed for this replication
set.seed(base_seed + rep_id)

# Helper functions
compute_G_in_true <- function(i, X0, Y0, Z0, tau) {
  n <- nrow(X0)
  p_cov <- nrow(Z0)
  d <- ncol(X0)

  s <- as.vector(X0[i, ] %*% t(Y0))
  w <- dpsi(s, tau = tau)

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

  G_net <- matrix(0, d, d)
  for (j in 1:n) {
    G_net <- G_net + w[j] * outer(Y_est[j, ], Y_est[j, ])
  }

  G_cov <- crossprod(Z_est)
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

get_precisions_ase_ose_1d <- function(i, X_target, eps_clip) {
  n <- length(X_target)
  idx_j <- setdiff(1:n, i)

  p_vals <- X_target * X_target[i]
  p_vals <- pmax(pmin(p_vals, 1 - eps_clip), eps_clip)

  # OSE Precision (Fisher Information)
  w_ose <- 1 / (p_vals[idx_j] * (1 - p_vals[idx_j]))
  prec_ose <- sum(X_target[idx_j]^2 * w_ose)

  # ASE Precision (Sandwich form)
  Delta <- sum(X_target^2)
  w_ase <- p_vals[idx_j] * (1 - p_vals[idx_j])
  M_mat <- sum(X_target[idx_j]^2 * w_ase)
  prec_ase <- Delta^2 / (M_mat + 1e-9)

  return(list(prec_ose = prec_ose, prec_ase = prec_ase))
}

ase_grdpg <- function(A, d) {
  spec <- eigen(A, symmetric = TRUE)
  X <- spec$vectors[, 1:d, drop = FALSE] %*% diag(sqrt(abs(spec$values[1:d])), nrow = d)
  return(list(X = X))
}

compute_ose_step_1d <- function(A, X_init, clipping_val) {
  n_nodes <- nrow(A)
  X_new <- matrix(0, n_nodes, 1)

  for (i in 1:n_nodes) {
    x_i <- X_init[i, 1]
    indices_j <- setdiff(1:n_nodes, i)

    p_i <- as.vector(X_init[indices_j, 1] * x_i)
    p_i <- pmax(pmin(p_i, 1 - clipping_val), clipping_val)

    resid <- A[i, indices_j] - p_i
    weight_score <- 1 / (p_i * (1 - p_i))
    grad <- sum(X_init[indices_j, 1] * (resid * weight_score))

    G <- sum(X_init[indices_j, 1]^2 * weight_score)

    step <- grad / (G + 1e-9)
    X_new[i, 1] <- x_i + step
  }
  return(X_new)
}

# Generate latent positions
cat(sprintf("Replication %d: Generating 1D latent positions...\n", rep_id))
i_vals <- 0:(n-1)
X0 <- matrix(0.75 * sin(pi * i_vals / (n - 1)) + 0.1, n, 1)

# For d=1, standard RDPG
S <- matrix(1, 1, 1)
Y0 <- X0 %*% S

# NO COVARIATE INFORMATION: Z0 = 0
Z0 <- matrix(0, p_cov, d)

P <- X0 %*% t(Y0)
cat(sprintf("Edge probability range: [%.4f, %.4f]\n\n", min(P), max(P)))

# Generate data
cat("Generating A and B...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

# B contains ONLY NOISE (no signal from Z0 %*% t(X0))
B <- matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# ===== METHOD 1: cgrdpg (Fisher scoring) =====
cat("Fitting cgrdpg model...\n")
cgrdpg_start <- Sys.time()
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (ncores <= 1) ncores <- max(1, parallel::detectCores() - 1)

fit <- tryCatch(
  fit_grdpg_cov_parallel(A, B, d = d, p = 1, q = 0,
                         maxit = maxit, tol = tol, tau = tau, ncores = ncores),
  error = function(e) fit_grdpg_cov(A, B, d = d, p = 1, q = 0,
                                     maxit = maxit, tol = tol, tau = tau)
)
cgrdpg_time <- as.numeric(difftime(Sys.time(), cgrdpg_start, units = "secs"))

# Procrustes for cgrdpg
X_cgrdpg_raw <- fit$X
M_cgrdpg <- t(X_cgrdpg_raw) %*% X0
svd_cgrdpg <- svd(M_cgrdpg)
Q_cgrdpg <- svd_cgrdpg$u %*% t(svd_cgrdpg$v)
X_cgrdpg <- X_cgrdpg_raw %*% Q_cgrdpg
Y_cgrdpg <- X_cgrdpg %*% S
Z_cgrdpg <- fit$Z

sse_cgrdpg <- sum((X_cgrdpg - X0)^2)
cat(sprintf("cgrdpg: converged=%s, iters=%d, time=%.1fs, SSE=%.4f\n",
            fit$converged, fit$iters, cgrdpg_time, sse_cgrdpg))

# ===== METHOD 2: ASE =====
cat("Computing ASE...\n")
ase_start <- Sys.time()
A_aug <- A
diag(A_aug) <- rowSums(A) / (n - 1)
ase_res <- ase_grdpg(A_aug, d = d)
X_ase_raw <- ase_res$X
ase_time <- as.numeric(difftime(Sys.time(), ase_start, units = "secs"))

# Procrustes for ASE
M_ase <- t(X_ase_raw) %*% X0
svd_ase <- svd(M_ase)
Q_ase <- svd_ase$u %*% t(svd_ase$v)
X_ase <- X_ase_raw %*% Q_ase

sse_ase <- sum((X_ase - X0)^2)
cat(sprintf("ASE: time=%.1fs, SSE=%.4f\n", ase_time, sse_ase))

# ===== METHOD 3: OSE =====
cat("Computing OSE...\n")
ose_step_start <- Sys.time()
X_ose_raw <- compute_ose_step_1d(A, X_ase_raw, clipping_val = eps_clip)
ose_step_time <- as.numeric(difftime(Sys.time(), ose_step_start, units = "secs"))
ose_time <- ase_time + ose_step_time

# Procrustes for OSE
M_ose <- t(X_ose_raw) %*% X0
svd_ose <- svd(M_ose)
Q_ose <- svd_ose$u %*% t(svd_ose$v)
X_ose <- X_ose_raw %*% Q_ose

sse_ose <- sum((X_ose - X0)^2)
cat(sprintf("OSE: time=%.1fs (ASE: %.1fs + step: %.1fs), SSE=%.4f\n\n",
            ose_time, ase_time, ose_step_time, sse_ose))

# ===== COVERAGE COMPUTATION =====
cat("Computing vertex-wise coverage...\n")
coverage_start <- Sys.time()
chi2_crit <- qchisq(0.95, df = d)

# cgrdpg coverage (TRUE and PLUGIN)
cgrdpg_true <- logical(n)
cgrdpg_plugin <- logical(n)

for (i in 1:n) {
  if (i %% 200 == 0) cat(sprintf("  Vertex %d/%d\n", i, n))

  # cgrdpg TRUE
  G_true <- compute_G_in_true(i, X0, Y0, Z0, tau)
  if (G_true[1,1] > 1e-10) {
    diff <- X0[i, 1] - X_cgrdpg[i, 1]
    mahal <- (n + p_cov) * diff^2 * G_true[1,1]
    cgrdpg_true[i] <- (mahal <= chi2_crit)
  } else {
    cgrdpg_true[i] <- NA
  }

  # cgrdpg PLUGIN
  G_plug <- compute_G_in_plugin(i, X_cgrdpg, Y_cgrdpg, Z_cgrdpg, tau)
  if (G_plug[1,1] > 1e-10) {
    mahal <- (n + p_cov) * diff^2 * G_plug[1,1]
    cgrdpg_plugin[i] <- (mahal <= chi2_crit)
  } else {
    cgrdpg_plugin[i] <- NA
  }
}

# ASE coverage (TRUE and PLUGIN)
ase_true <- logical(n)
ase_plugin <- logical(n)

for (i in 1:n) {
  prec_true <- get_precisions_ase_ose_1d(i, X0[,1], eps_clip)
  prec_plug <- get_precisions_ase_ose_1d(i, X_ase[,1], eps_clip)

  diff <- X0[i, 1] - X_ase[i, 1]

  if (prec_true$prec_ase > 1e-10) {
    mahal_true <- diff^2 * prec_true$prec_ase
    ase_true[i] <- (mahal_true <= chi2_crit)
  } else {
    ase_true[i] <- NA
  }

  if (prec_plug$prec_ase > 1e-10) {
    mahal_plug <- diff^2 * prec_plug$prec_ase
    ase_plugin[i] <- (mahal_plug <= chi2_crit)
  } else {
    ase_plugin[i] <- NA
  }
}

# OSE coverage (TRUE and PLUGIN)
ose_true <- logical(n)
ose_plugin <- logical(n)

for (i in 1:n) {
  prec_true <- get_precisions_ase_ose_1d(i, X0[,1], eps_clip)
  prec_plug <- get_precisions_ase_ose_1d(i, X_ose[,1], eps_clip)

  diff <- X0[i, 1] - X_ose[i, 1]

  if (prec_true$prec_ose > 1e-10) {
    mahal_true <- diff^2 * prec_true$prec_ose
    ose_true[i] <- (mahal_true <= chi2_crit)
  } else {
    ose_true[i] <- NA
  }

  if (prec_plug$prec_ose > 1e-10) {
    mahal_plug <- diff^2 * prec_plug$prec_ose
    ose_plugin[i] <- (mahal_plug <= chi2_crit)
  } else {
    ose_plugin[i] <- NA
  }
}

coverage_time <- as.numeric(difftime(Sys.time(), coverage_start, units = "secs"))

# Compile results
result <- list(
  rep_id = rep_id,
  seed = base_seed + rep_id,
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
  S = S,
  sse = list(
    cgrdpg = sse_cgrdpg,
    ase = sse_ase,
    ose = sse_ose
  ),
  coverage = list(
    cgrdpg_true = cgrdpg_true,
    cgrdpg_plugin = cgrdpg_plugin,
    ase_true = ase_true,
    ase_plugin = ase_plugin,
    ose_true = ose_true,
    ose_plugin = ose_plugin
  ),
  overall_cov = c(
    cgrdpg_true = mean(cgrdpg_true, na.rm = TRUE),
    cgrdpg_plugin = mean(cgrdpg_plugin, na.rm = TRUE),
    ase_true = mean(ase_true, na.rm = TRUE),
    ase_plugin = mean(ase_plugin, na.rm = TRUE),
    ose_true = mean(ose_true, na.rm = TRUE),
    ose_plugin = mean(ose_plugin, na.rm = TRUE)
  ),
  n_na = c(
    cgrdpg_true = sum(is.na(cgrdpg_true)),
    cgrdpg_plugin = sum(is.na(cgrdpg_plugin)),
    ase_true = sum(is.na(ase_true)),
    ase_plugin = sum(is.na(ase_plugin)),
    ose_true = sum(is.na(ose_true)),
    ose_plugin = sum(is.na(ose_plugin))
  ),
  timing = list(
    cgrdpg_time = cgrdpg_time,
    ase_time = ase_time,
    ose_time = ose_time,
    coverage_time = coverage_time,
    rep_time_min = as.numeric(difftime(Sys.time(), cgrdpg_start, units = "mins"))
  ),
  converged = fit$converged,
  iterations = fit$iters
)

# Save result
output_dir <- "results_1d_no_info_n1000"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(output_dir, sprintf("rep_%03d.rds", rep_id))
saveRDS(result, output_file)

cat("\n============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETE\n", rep_id))
cat("============================================================================\n")
cat(sprintf("SSE:     cgrdpg=%.4f, ASE=%.4f, OSE=%.4f\n", sse_cgrdpg, sse_ase, sse_ose))
cat(sprintf("Coverage (TRUE):   cgrdpg=%.1f%%, ASE=%.1f%%, OSE=%.1f%%\n",
            100*mean(cgrdpg_true, na.rm=TRUE),
            100*mean(ase_true, na.rm=TRUE),
            100*mean(ose_true, na.rm=TRUE)))
cat(sprintf("Coverage (PLUGIN): cgrdpg=%.1f%%, ASE=%.1f%%, OSE=%.1f%%\n",
            100*mean(cgrdpg_plugin, na.rm=TRUE),
            100*mean(ase_plugin, na.rm=TRUE),
            100*mean(ose_plugin, na.rm=TRUE)))
cat(sprintf("Timing: cgrdpg=%.1fs, ASE=%.1fs, OSE=%.1fs, Coverage=%.1fs\n",
            cgrdpg_time, ase_time, ose_time, coverage_time))
cat(sprintf("Result saved to: %s\n", output_file))
cat("============================================================================\n")
