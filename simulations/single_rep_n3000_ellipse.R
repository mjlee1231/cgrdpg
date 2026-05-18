# Single replication script for n=3000 ellipse coverage
# Uses SLURM_ARRAY_TASK_ID as the replication number and seed
library(cgrdpg)

# Get replication number from SLURM array task ID
rep_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

generate_latent_curve <- function(n) {
  t_vals <- (1:n) / n
  X <- matrix(0, n, 3)
  X[, 1] <- 0.12 * sin(2 * pi * t_vals) + 0.55
  X[, 2] <- 0.12 * cos(2 * pi * t_vals) + 0.55
  X[, 3] <- 0.08 * cos(4 * pi * t_vals)
  return(X)
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

  # G_in with normalization as per theoretical formula
  G_in <- (G_net + G_cov) / (n + p_cov)

  return(G_in)
}

cat("============================================================================\n")
cat(sprintf("  Single Replication: n=3000, Rep ID=%d\n", rep_id))
cat("============================================================================\n\n")

# Set parameters
set.seed(1000 + rep_id)  # Different seed for each replication
n <- 3000
p_cov <- 1500  # p_cov = n/2
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

# Start timing
rep_start <- Sys.time()

# Generate data
cat("1. Generating adjacency matrix A...\n")
A <- (runif(n = n^2, min = 0, max = 1) < P) * 1.0
A <- A * upper.tri(x = A, diag = FALSE) + t(A * upper.tri(x = A, diag = FALSE))

cat("2. Generating covariate matrix B...\n")
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)

# Fit model
cat("3. Fitting GRDPG model...\n")
fit_start <- Sys.time()
fit <- fit_grdpg_cov(A, B, d = d, p = 2, q = 1,
                      maxit = maxit, tol = tol, tau = tau)
fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
cat(sprintf("   Model fitting: %.2f seconds (%.2f minutes)\n", fit_time, fit_time/60))

# Procrustes alignment
cat("4. Procrustes alignment...\n")
X_est <- fit$X
Z_est <- fit$Z
M <- t(X_est) %*% X0
svd_res <- svd(M)
R <- svd_res$u %*% t(svd_res$v)
X_aligned <- X_est %*% R
Z_updated <- B %*% X_aligned %*% solve(t(X_aligned) %*% X_aligned)

SSE <- sum((X_aligned - X0)^2)

# Check ALL nodes for ellipse coverage using PLUG-IN G_in
cat(sprintf("5. Computing coverage for all %d nodes (plug-in G_in)...\n", n))
cov_start <- Sys.time()

chi2_crit <- qchisq(0.95, df = d)
Y_aligned <- X_aligned %*% S
in_ellipse <- logical(n)

for (i in 1:n) {
  if (i %% 500 == 0) cat(sprintf("  Node %d/%d\n", i, n))

  # Compute PLUG-IN G_in using estimated parameters
  G_in_plugin <- compute_G_in_plugin(i, X_aligned, Y_aligned, Z_updated, tau = tau)

  eig_vals <- eigen(G_in_plugin, only.values = TRUE)$values
  if (min(eig_vals) < 1e-10) {
    in_ellipse[i] <- NA
    next
  }

  diff <- X0[i, ] - X_aligned[i, ]

  # Mahalanobis distance (scaled by n+p_cov)
  mahal_dist <- (n + p_cov) * as.numeric(t(diff) %*% G_in_plugin %*% diff)

  in_ellipse[i] <- (mahal_dist <= chi2_crit)
}

cov_time <- as.numeric(difftime(Sys.time(), cov_start, units = "secs"))
rep_time <- as.numeric(difftime(Sys.time(), rep_start, units = "secs"))

coverage <- mean(in_ellipse, na.rm = TRUE)
n_inside <- sum(in_ellipse, na.rm = TRUE)
n_outside <- sum(!in_ellipse, na.rm = TRUE)
n_na <- sum(is.na(in_ellipse))

cat("\n============================================================================\n")
cat(sprintf("  REPLICATION %d COMPLETE\n", rep_id))
cat("============================================================================\n\n")

cat("--- Results ---\n")
cat(sprintf("Coverage:           %.1f%% (%d/%d nodes)\n", 100*coverage, n_inside, n))
cat(sprintf("Outside ellipse:    %d nodes\n", n_outside))
cat(sprintf("NA (numerical):     %d nodes\n", n_na))
cat(sprintf("SSE:                %.4f\n", SSE))
cat(sprintf("Total time:         %.2f seconds (%.2f minutes) (%.2f hours)\n",
            rep_time, rep_time/60, rep_time/3600))
cat(sprintf("  - Model fitting:  %.2f seconds (%.2f minutes)\n", fit_time, fit_time/60))
cat(sprintf("  - Coverage check: %.2f seconds\n\n", cov_time))

# Create output directory if it doesn't exist
output_dir <- "results_n3000"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save results to unique file in output directory
output_file <- file.path(output_dir, sprintf("results_n3000_rep%03d.rds", rep_id))
results <- list(
  rep_id = rep_id,
  coverage = coverage,
  n_inside = n_inside,
  n_outside = n_outside,
  n_na = n_na,
  SSE = SSE,
  rep_time = rep_time,
  fit_time = fit_time,
  cov_time = cov_time,
  tau = tau,
  min_p = min_p,
  max_p = max_p,
  n = n,
  p_cov = p_cov
)

saveRDS(results, output_file)
cat(sprintf("Results saved to: %s\n", output_file))
cat("============================================================================\n")
