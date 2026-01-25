# Diagnostic: Compare TRUE G_in vs PLUG-IN G_in coverage at n=500
# This will show if plug-in estimation is the bottleneck
library(cgrdpg)

generate_latent_curve <- function(n) {
  t_vals <- (1:n) / n
  X <- matrix(0, n, 3)
  X[, 1] <- 0.12 * sin(2 * pi * t_vals) + 0.55
  X[, 2] <- 0.12 * cos(2 * pi * t_vals) + 0.55
  X[, 3] <- 0.08 * cos(4 * pi * t_vals)
  return(X)
}

compute_G_in_true <- function(i, X0, Y0, Z0, tau) {
  # Uses TRUE parameters (oracle)
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

  # CORRECTED: Include normalization
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

compute_G_in_plugin <- function(i, X_est, Y_est, Z_est, tau) {
  # Uses ESTIMATED parameters (plug-in)
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

  # CORRECTED: Include normalization
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

cat("============================================================================\n")
cat("  DIAGNOSTIC: TRUE G_in vs PLUG-IN G_in at n=500\n")
cat("============================================================================\n\n")

set.seed(598)
n <- 500
p_cov <- 250
d <- 3
maxit <- 3
tol <- 0.01

# Generate true parameters
cat("Generating true parameters...\n")
X0 <- generate_latent_curve(n)
S <- diag(c(1, 1, -1))
Y0 <- X0 %*% S
Z0 <- matrix(rnorm(p_cov * d), p_cov, d)

# Compute edge probability matrix
P <- X0 %*% t(Y0)
min_p <- min(P)
max_p <- max(P)
buffer <- 0.0001
tau_auto <- min(min_p, 1 - max_p) - buffer
tau <- min(1e-3, tau_auto)

cat(sprintf("Edge Probability Range: [%.4f, %.4f]\n", min_p, max_p))
cat(sprintf("Selected tau: %.6f\n\n", tau))

# Generate data
cat("Generating adjacency matrix A...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

cat("Generating covariate matrix B...\n")
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# Fit model
cat("Fitting GRDPG model...\n")
fit_start <- Sys.time()
fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 1,
                      maxit = maxit, tol = tol, tau = tau)
fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
cat(sprintf("Model fitting: %.2f seconds (%.2f minutes)\n\n", fit_time, fit_time/60))

# Procrustes alignment
cat("Procrustes alignment...\n")
X_est <- fit$X
Z_est <- fit$Z

Uhat <- svd(X_est)$u
U0 <- svd(X0)$u
W_svd <- svd(t(Uhat) %*% U0)
R <- W_svd$u %*% t(W_svd$v)
X_aligned <- X_est %*% R
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)


SSE <- sum((X_aligned - X0)^2)
cat(sprintf("SSE: %.4f\n\n", SSE))

# Check coverage with BOTH methods
cat("Computing coverage with TRUE G_in (oracle)...\n")
chi2_crit <- qchisq(0.95, df = d)
Y_aligned <- X_aligned %*% S
in_ellipse_true <- logical(n)

for (i in 1:n) {
  if (i %% 100 == 0) cat(sprintf("  Node %d/%d\n", i, n))

  # Use TRUE G_in
  G_in_true <- compute_G_in_true(i, X0, Y0, Z0, tau = tau)

  eig_vals <- eigen(G_in_true, only.values = TRUE)$values
  if (min(eig_vals) < 1e-10) {
    in_ellipse_true[i] <- NA
    next
  }

  diff <- X0[i, ] - X_aligned[i, ]
  mahal_dist <- (n + p_cov) * as.numeric(t(diff) %*% G_in_true %*% diff)

  in_ellipse_true[i] <- (mahal_dist <= chi2_crit)
}

cat("\nComputing coverage with PLUG-IN G_in (estimated)...\n")
in_ellipse_plugin <- logical(n)

for (i in 1:n) {
  if (i %% 100 == 0) cat(sprintf("  Node %d/%d\n", i, n))

  # Use PLUG-IN G_in
  G_in_plugin <- compute_G_in_plugin(i, X_aligned, Y_aligned, Z_updated, tau = tau)

  eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
  if (min(eig_vals) < 1e-10) {
    in_ellipse_plugin[i] <- NA
    next
  }

  diff <- X0[i, ] - X_aligned[i, ]
  mahal_dist <- (n + p_cov) * as.numeric(t(diff) %*% G_in_plugin %*% diff)

  in_ellipse_plugin[i] <- (mahal_dist <= chi2_crit)
}

# Results
coverage_true <- mean(in_ellipse_true, na.rm = TRUE)
coverage_plugin <- mean(in_ellipse_plugin, na.rm = TRUE)
n_inside_true <- sum(in_ellipse_true, na.rm = TRUE)
n_inside_plugin <- sum(in_ellipse_plugin, na.rm = TRUE)
n_na_true <- sum(is.na(in_ellipse_true))
n_na_plugin <- sum(is.na(in_ellipse_plugin))

cat("\n============================================================================\n")
cat("  COMPARISON RESULTS\n")
cat("============================================================================\n\n")

cat(sprintf("n = %d, p_cov = %d, d = %d\n", n, p_cov, d))
cat(sprintf("SSE = %.4f\n\n", SSE))

cat("--- TRUE G_in (Oracle with X0, Y0, Z0) ---\n")
cat(sprintf("Coverage:        %.1f%% (%d/%d nodes)\n", 100*coverage_true, n_inside_true, n))
cat(sprintf("Outside:         %d nodes\n", sum(!in_ellipse_true, na.rm = TRUE)))
cat(sprintf("NAs:             %d nodes\n\n", n_na_true))

cat("--- PLUG-IN G_in (Estimated with X_aligned, Y_aligned, Z_updated) ---\n")
cat(sprintf("Coverage:        %.1f%% (%d/%d nodes)\n", 100*coverage_plugin, n_inside_plugin, n))
cat(sprintf("Outside:         %d nodes\n", sum(!in_ellipse_plugin, na.rm = TRUE)))
cat(sprintf("NAs:             %d nodes\n\n", n_na_plugin))

cat("--- Analysis ---\n")
cat(sprintf("Coverage difference: %.1f%% (TRUE) - %.1f%% (PLUG-IN) = %.1f%%\n",
            100*coverage_true, 100*coverage_plugin, 100*(coverage_true - coverage_plugin)))

if (coverage_true >= 0.90 && coverage_plugin < 0.70) {
  cat("\n*** INTERPRETATION ***\n")
  cat("TRUE G_in gives ~95%% coverage (as expected by theory)\n")
  cat("PLUG-IN G_in gives low coverage\n")
  cat("=> The bottleneck is PLUG-IN ESTIMATION, not the theory itself\n")
  cat("=> Need oracle rotation Q or better variance estimation\n")
} else if (coverage_true < 0.70 && coverage_plugin < 0.70) {
  cat("\n*** INTERPRETATION ***\n")
  cat("BOTH TRUE and PLUG-IN give low coverage\n")
  cat("=> Problem is deeper than plug-in estimation\n")
  cat("=> Check: Procrustes rotation, theory assumptions, or implementation\n")
} else {
  cat("\n*** INTERPRETATION ***\n")
  cat("Results suggest moderate plug-in estimation penalty\n")
}

cat("============================================================================\n")
