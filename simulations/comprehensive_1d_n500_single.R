#!/usr/bin/env Rscript
# Comprehensive 1D simulation: Fisher vs ASE vs OSE (SINGLE REPLICATION)
# Latent position: X_i = 0.75 * sin(π * i/(n-1)) + 0.1
# Designed for SLURM array: each job runs one replication

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript comprehensive_1d_n500_single.R <rep_id>")
}

rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 100) {
  stop("rep_id must be between 1 and 100")
}

library(cgrdpg)

cat("============================================================================\n")
cat(sprintf("  1D SIMULATION - REPLICATION %d/100\n", rep_id))
cat("  Comparing: FISHER (with covariates) vs ASE vs OSE (no covariates)\n")
cat("============================================================================\n\n")

# Fixed parameters
n <- 500
p_cov <- 250
d <- 1
maxit <- 3
tol <- 0.01
base_seed <- 598
eps_clip <- 1e-3

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
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

# Compute tau
P <- X0 %*% t(Y0)
min_p <- min(P)
max_p <- max(P)
buffer <- 0.0001
tau_auto <- min(min_p, 1 - max_p) - buffer
tau <- min(1e-3, tau_auto)

# Generate data
cat("Generating A and B...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# ===== METHOD 1: FISHER-SCORING =====
cat("Fitting FISHER model...\n")
fisher_start <- Sys.time()
fit <- fit_grdpg_cov(A, B, d = d, p = 1, q = 0,
                      maxit = maxit, tol = tol, tau = tau)
fisher_time <- as.numeric(difftime(Sys.time(), fisher_start, units = "secs"))

# Procrustes for Fisher
X_fisher_raw <- fit$X
M_fisher <- t(X_fisher_raw) %*% X0
svd_fisher <- svd(M_fisher)
Q_fisher <- svd_fisher$u %*% t(svd_fisher$v)
X_fisher <- X_fisher_raw %*% Q_fisher
Y_fisher <- X_fisher %*% S
Z_fisher <- fit$Z

sse_fisher <- sum((X_fisher - X0)^2)

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

# ===== METHOD 3: OSE =====
cat("Computing OSE...\n")
ose_start <- Sys.time()
X_ose_raw <- compute_ose_step_1d(A, X_ase_raw, clipping_val = eps_clip)
ose_time <- as.numeric(difftime(Sys.time(), ose_start, units = "secs"))

# Procrustes for OSE
M_ose <- t(X_ose_raw) %*% X0
svd_ose <- svd(M_ose)
Q_ose <- svd_ose$u %*% t(svd_ose$v)
X_ose <- X_ose_raw %*% Q_ose

sse_ose <- sum((X_ose - X0)^2)

# ===== COVERAGE COMPUTATION =====
cat("Computing coverage...\n")
coverage_start <- Sys.time()
chi2_crit <- qchisq(0.95, df = d)

# Fisher coverage (TRUE and PLUGIN)
fisher_true <- logical(n)
fisher_plugin <- logical(n)

for (i in 1:n) {
  # Fisher TRUE
  G_true <- compute_G_in_true(i, X0, Y0, Z0, tau)
  if (G_true[1,1] > 1e-10) {
    diff <- X0[i, 1] - X_fisher[i, 1]
    mahal <- (n + p_cov) * diff^2 * G_true[1,1]
    fisher_true[i] <- (mahal <= chi2_crit)
  } else {
    fisher_true[i] <- NA
  }

  # Fisher PLUGIN
  G_plug <- compute_G_in_plugin(i, X_fisher, Y_fisher, Z_fisher, tau)
  if (G_plug[1,1] > 1e-10) {
    mahal <- (n + p_cov) * diff^2 * G_plug[1,1]
    fisher_plugin[i] <- (mahal <= chi2_crit)
  } else {
    fisher_plugin[i] <- NA
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
  n = n,
  p_cov = p_cov,
  d = d,
  sse = list(
    fisher = sse_fisher,
    ase = sse_ase,
    ose = sse_ose
  ),
  coverage = list(
    fisher_true = fisher_true,
    fisher_plugin = fisher_plugin,
    ase_true = ase_true,
    ase_plugin = ase_plugin,
    ose_true = ose_true,
    ose_plugin = ose_plugin
  ),
  timing = list(
    fisher = fisher_time,
    ase = ase_time,
    ose = ose_time,
    coverage = coverage_time
  )
)

# Save result
output_dir <- "comprehensive_1d_n500_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

output_file <- file.path(output_dir, sprintf("rep_%03d.rds", rep_id))
saveRDS(result, output_file)

cat("\n============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETE\n", rep_id))
cat("============================================================================\n")
cat(sprintf("SSE:     Fisher=%.4f, ASE=%.4f, OSE=%.4f\n", sse_fisher, sse_ase, sse_ose))
cat(sprintf("Coverage (TRUE):   Fisher=%.1f%%, ASE=%.1f%%, OSE=%.1f%%\n",
            100*mean(fisher_true, na.rm=TRUE),
            100*mean(ase_true, na.rm=TRUE),
            100*mean(ose_true, na.rm=TRUE)))
cat(sprintf("Coverage (PLUGIN): Fisher=%.1f%%, ASE=%.1f%%, OSE=%.1f%%\n",
            100*mean(fisher_plugin, na.rm=TRUE),
            100*mean(ase_plugin, na.rm=TRUE),
            100*mean(ose_plugin, na.rm=TRUE)))
cat(sprintf("Timing: Fisher=%.1fs, ASE=%.1fs, OSE=%.1fs, Coverage=%.1fs\n",
            fisher_time, ase_time, ose_time, coverage_time))
cat(sprintf("Result saved to: %s\n", output_file))
cat("============================================================================\n")
