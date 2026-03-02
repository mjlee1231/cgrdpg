#!/usr/bin/env Rscript
# Comprehensive 1D simulation: Fisher vs ASE vs OSE (100 replications)
# Latent position: X_i = 0.75 * sin(π * i/(n-1)) + 0.1
# Compares all three methods on the SAME datasets

library(cgrdpg)

cat("============================================================================\n")
cat("  COMPREHENSIVE 1D SIMULATION (n=500, 100 replications)\n")
cat("  Comparing: FISHER (with covariates) vs ASE vs OSE (no covariates)\n")
cat("============================================================================\n\n")

# Fixed parameters
n <- 500
p_cov <- 250
d <- 1
n_reps <- 100
maxit <- 30
tol <- 0.01
base_seed <- 598

# Latent position specification (per user request)
cat("Latent position: X_i = 0.75 * sin(π * i/(n-1)) + 0.1\n\n")

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

# Storage for results
results <- vector("list", n_reps)

# Progress tracking
global_start <- Sys.time()
cat(sprintf("Starting %d replications...\n\n", n_reps))

for (rep in 1:n_reps) {
  rep_start <- Sys.time()
  set.seed(base_seed + rep)

  if (rep %% 10 == 0) cat(sprintf("Replication %d/%d\n", rep, n_reps))

  # Generate 1D latent positions: 0.75 * sin(π * i/(n-1)) + 0.1
  i_vals <- 1:n
  X0 <- matrix(0.75 * sin(pi * i_vals / (n - 1)) + 0.1, n, 1)

  # Positive signature for d=1
  S <- matrix(1, 1, 1)
  Y0 <- X0 %*% S
  Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

  # Edge probabilities
  P <- X0 %*% t(Y0)
  min_p <- min(P)
  max_p <- max(P)
  buffer <- 0.0001
  tau_auto <- min(min_p, 1 - max_p) - buffer
  tau <- min(1e-3, tau_auto)
  eps_clip <- 1e-10

  # Generate data (same for all methods)
  A <- (runif(n^2) < P) * 1.0
  A <- A * upper.tri(A, diag = FALSE) + t(A * upper.tri(A, diag = FALSE))
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

  # ===== METHOD 1: FISHER-SCORING (with covariates) =====
  fisher_start <- Sys.time()
  fit <- fit_grdpg_cov(A, B, d = d, p = 1, q = 0,
                       maxit = maxit, tol = tol, tau = tau)
  fisher_time <- as.numeric(difftime(Sys.time(), fisher_start, units = "secs"))

  # Procrustes alignment
  M <- t(fit$X) %*% X0
  svd_res <- svd(M)
  R <- svd_res$u %*% t(svd_res$v)
  X_fisher <- fit$X %*% R
  Z_fisher <- B %*% X_fisher %*% solve(t(X_fisher) %*% X_fisher)
  Y_fisher <- X_fisher %*% S

  sse_fisher <- sum((X_fisher - X0)^2)

  # ===== METHOD 2: ASE (no covariates) =====
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

  # ===== METHOD 3: OSE (no covariates, 1-step from ASE) =====
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
  rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "secs"))

  # Store results
  results[[rep]] <- list(
    rep_id = rep,
    seed = base_seed + rep,
    sse = list(fisher = sse_fisher, ase = sse_ase, ose = sse_ose),
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
      coverage = coverage_time,
      total = rep_time
    )
  )
}

global_time <- as.numeric(difftime(Sys.time(), global_start, units = "secs"))

cat(sprintf("\n\nAll %d replications completed in %.1f seconds (%.1f minutes)\n\n",
            n_reps, global_time, global_time/60))

# ===== AGGREGATE RESULTS =====
cat("============================================================================\n")
cat("  AGGREGATING RESULTS ACROSS 100 REPLICATIONS\n")
cat("============================================================================\n\n")

# Extract coverage matrices
fisher_true_mat <- matrix(NA, n, n_reps)
fisher_plugin_mat <- matrix(NA, n, n_reps)
ase_true_mat <- matrix(NA, n, n_reps)
ase_plugin_mat <- matrix(NA, n, n_reps)
ose_true_mat <- matrix(NA, n, n_reps)
ose_plugin_mat <- matrix(NA, n, n_reps)

sse_fisher <- numeric(n_reps)
sse_ase <- numeric(n_reps)
sse_ose <- numeric(n_reps)

time_fisher <- numeric(n_reps)
time_ase <- numeric(n_reps)
time_ose <- numeric(n_reps)
time_coverage <- numeric(n_reps)
time_total <- numeric(n_reps)

for (r in 1:n_reps) {
  fisher_true_mat[, r] <- results[[r]]$coverage$fisher_true
  fisher_plugin_mat[, r] <- results[[r]]$coverage$fisher_plugin
  ase_true_mat[, r] <- results[[r]]$coverage$ase_true
  ase_plugin_mat[, r] <- results[[r]]$coverage$ase_plugin
  ose_true_mat[, r] <- results[[r]]$coverage$ose_true
  ose_plugin_mat[, r] <- results[[r]]$coverage$ose_plugin

  sse_fisher[r] <- results[[r]]$sse$fisher
  sse_ase[r] <- results[[r]]$sse$ase
  sse_ose[r] <- results[[r]]$sse$ose

  time_fisher[r] <- results[[r]]$timing$fisher
  time_ase[r] <- results[[r]]$timing$ase
  time_ose[r] <- results[[r]]$timing$ose
  time_coverage[r] <- results[[r]]$timing$coverage
  time_total[r] <- results[[r]]$timing$total
}

# Compute vertex-wise coverage
vertex_cov_fisher_true <- rowMeans(fisher_true_mat, na.rm = TRUE)
vertex_cov_fisher_plugin <- rowMeans(fisher_plugin_mat, na.rm = TRUE)
vertex_cov_ase_true <- rowMeans(ase_true_mat, na.rm = TRUE)
vertex_cov_ase_plugin <- rowMeans(ase_plugin_mat, na.rm = TRUE)
vertex_cov_ose_true <- rowMeans(ose_true_mat, na.rm = TRUE)
vertex_cov_ose_plugin <- rowMeans(ose_plugin_mat, na.rm = TRUE)

# Print summary
cat("VERTEX-WISE COVERAGE (mean across vertices and replications):\n\n")
cat(sprintf("FISHER-TRUE:   %.2f%%\n", 100*mean(vertex_cov_fisher_true, na.rm=TRUE)))
cat(sprintf("FISHER-PLUGIN: %.2f%%\n", 100*mean(vertex_cov_fisher_plugin, na.rm=TRUE)))
cat(sprintf("ASE-TRUE:      %.2f%%\n", 100*mean(vertex_cov_ase_true, na.rm=TRUE)))
cat(sprintf("ASE-PLUGIN:    %.2f%%\n", 100*mean(vertex_cov_ase_plugin, na.rm=TRUE)))
cat(sprintf("OSE-TRUE:      %.2f%%\n", 100*mean(vertex_cov_ose_true, na.rm=TRUE)))
cat(sprintf("OSE-PLUGIN:    %.2f%%\n\n", 100*mean(vertex_cov_ose_plugin, na.rm=TRUE)))

cat("SSE (mean across replications):\n\n")
cat(sprintf("FISHER: %.4f (SD: %.4f)\n", mean(sse_fisher), sd(sse_fisher)))
cat(sprintf("ASE:    %.4f (SD: %.4f)\n", mean(sse_ase), sd(sse_ase)))
cat(sprintf("OSE:    %.4f (SD: %.4f)\n\n", mean(sse_ose), sd(sse_ose)))

cat("TIMING (mean per replication in seconds):\n\n")
cat(sprintf("FISHER:         %.2f sec\n", mean(time_fisher)))
cat(sprintf("ASE:            %.2f sec\n", mean(time_ase)))
cat(sprintf("OSE:            %.2f sec\n", mean(time_ose)))
cat(sprintf("Coverage comp:  %.2f sec\n", mean(time_coverage)))
cat(sprintf("Total per rep:  %.2f sec\n\n", mean(time_total)))

# Save results
output <- list(
  n = n,
  p_cov = p_cov,
  d = d,
  n_reps = n_reps,
  vertex_coverage = list(
    fisher_true = vertex_cov_fisher_true,
    fisher_plugin = vertex_cov_fisher_plugin,
    ase_true = vertex_cov_ase_true,
    ase_plugin = vertex_cov_ase_plugin,
    ose_true = vertex_cov_ose_true,
    ose_plugin = vertex_cov_ose_plugin
  ),
  sse = list(
    fisher = sse_fisher,
    ase = sse_ase,
    ose = sse_ose
  ),
  timing = list(
    fisher = time_fisher,
    ase = time_ase,
    ose = time_ose,
    coverage = time_coverage,
    total = time_total
  ),
  all_results = results
)

saveRDS(output, "comprehensive_1d_n500_100reps_results.rds")
cat("Results saved to: comprehensive_1d_n500_100reps_results.rds\n")

cat("\n============================================================================\n")
cat("  SIMULATION COMPLETE\n")
cat("============================================================================\n")
