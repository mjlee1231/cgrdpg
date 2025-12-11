# Visualization of Asymptotic Normality
# Create scatter plots showing the distribution of estimated positions
# across repeated experiments for randomly selected latent positions

library(cgrdpg)
library(ggplot2)
library(reshape2)

# ------------------------------------------------------------------------------
# Data Generation Functions
# ------------------------------------------------------------------------------
generate_latent_curve <- function(n) {
  t_vals <- (1:n) / n
  X <- matrix(0, n, 3)
  X[, 1] <- 0.15 * sin(2 * pi * t_vals) + 0.6
  X[, 2] <- 0.15 * cos(2 * pi * t_vals) + 0.6
  X[, 3] <- 0.15 * cos(4 * pi * t_vals)
  return(X)
}

generate_grdpg_data <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1) {
  X0 <- generate_latent_curve(n)
  S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))
  Y0 <- X0 %*% S
  P <- X0 %*% t(Y0)

  A <- matrix(rbinom(n * n, 1, as.vector(P)), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0

  Z0 <- matrix(rnorm(p_cov * d, sd = 0.5), p_cov, d)
  B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 0.3), p_cov, n)

  list(A = A, B = B, X0 = X0, Y0 = Y0, Z0 = Z0, P = P, S = S)
}

compute_G_in <- function(i, X0, Y0, Z0, tau = 0.05) {
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

procrustes_rotation <- function(X_hat, X0) {
  M <- t(X_hat) %*% X0
  svd_res <- svd(M)
  R <- svd_res$u %*% t(svd_res$v)
  return(R)
}

# ------------------------------------------------------------------------------
# Collect Estimates Across Replicates
# ------------------------------------------------------------------------------
collect_estimates <- function(n, p_cov, d = 3, p_sig = 2, q_sig = 1,
                               n_reps = 500, vertices_to_track = NULL,
                               seed = 456, verbose = TRUE) {
  set.seed(seed)

  # Start timing
  start_time <- Sys.time()

  # Generate one dataset to get true positions
  data_ref <- generate_grdpg_data(n, p_cov, d, p_sig, q_sig)
  X0_ref <- data_ref$X0
  Y0_ref <- data_ref$Y0
  Z0_ref <- data_ref$Z0

  if (is.null(vertices_to_track)) {
    # Track a few vertices across the curve
    vertices_to_track <- c(5, 15, 25, 35, 45)
  }

  n_vertices <- length(vertices_to_track)

  # Storage: for each vertex, store estimates across reps
  estimates <- array(0, dim = c(n_reps, n_vertices, d))
  dimnames(estimates)[[2]] <- paste0("vertex_", vertices_to_track)
  dimnames(estimates)[[3]] <- paste0("coord_", 1:d)

  # Also compute theoretical covariances
  theoretical_cov <- array(0, dim = c(n_vertices, d, d))

  if (verbose) {
    cat("=== Collecting Estimates for Asymptotic Normality ===\n")
    cat(sprintf("n=%d, p_cov=%d, reps=%d\n", n, p_cov, n_reps))
    cat(sprintf("Tracking %d vertices: %s\n",
                n_vertices, paste(vertices_to_track, collapse = ", ")))
    cat(sprintf("Start time: %s\n\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  }

  replicate_times <- numeric(n_reps)

  for (rep in 1:n_reps) {
    rep_start <- Sys.time()

    if (verbose && rep %% 50 == 0) {
      elapsed <- as.numeric(difftime(rep_start, start_time, units = "secs"))
      eta <- (elapsed / rep) * (n_reps - rep)
      cat(sprintf("Replicate %d/%d (Elapsed: %.1fs, ETA: %.1fs)\n",
                  rep, n_reps, elapsed, eta))
    }

    # Generate new data with SAME X0 structure
    data <- generate_grdpg_data(n, p_cov, d, p_sig, q_sig)

    # Fit model
    fit <- fit_grdpg_cov(data$A, data$B, d = d, p = p_sig, q = q_sig,
                         maxit = 10, tol = 1e-3)

    # Align to reference X0 (always the same curve)
    R <- procrustes_rotation(fit$X, X0_ref)
    X_aligned <- fit$X %*% R

    # Store estimates
    for (idx in 1:n_vertices) {
      i <- vertices_to_track[idx]
      estimates[rep, idx, ] <- X_aligned[i, ]
    }

    # Record time for this replicate
    replicate_times[rep] <- as.numeric(difftime(Sys.time(), rep_start, units = "secs"))
  }

  # End timing
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Compute theoretical covariance for each tracked vertex
  for (idx in 1:n_vertices) {
    i <- vertices_to_track[idx]
    G_in <- compute_G_in(i, X0_ref, Y0_ref, Z0_ref, tau = 0.05)
    theoretical_cov[idx, , ] <- solve(G_in)
  }

  if (verbose) {
    cat("\n=== Collection Complete ===\n")
    cat(sprintf("Total time: %.2f seconds (%.2f minutes)\n",
                total_time, total_time / 60))
    cat(sprintf("Average time per replicate: %.3f seconds\n",
                mean(replicate_times)))
    cat(sprintf("End time: %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
  }

  list(
    estimates = estimates,
    X0_ref = X0_ref,
    vertices_tracked = vertices_to_track,
    theoretical_cov = theoretical_cov,
    n_reps = n_reps,
    timing = list(
      total_time_seconds = total_time,
      mean_replicate_time = mean(replicate_times),
      median_replicate_time = median(replicate_times),
      all_replicate_times = replicate_times,
      start_time = start_time,
      end_time = end_time
    )
  )
}

# ------------------------------------------------------------------------------
# Visualization Functions
# ------------------------------------------------------------------------------
plot_scatter_2d <- function(collected_data, vertex_idx = 1, coords = c(1, 2)) {
  vertices <- collected_data$vertices_tracked
  vertex_id <- vertices[vertex_idx]
  X0 <- collected_data$X0_ref

  # Extract estimates
  estimates <- collected_data$estimates[, vertex_idx, ]
  x_est <- estimates[, coords[1]]
  y_est <- estimates[, coords[2]]

  # True value
  x_true <- X0[vertex_id, coords[1]]
  y_true <- X0[vertex_id, coords[2]]

  # Theoretical covariance
  Sigma <- collected_data$theoretical_cov[vertex_idx, , ]
  Sigma_2d <- Sigma[coords, coords]

  # Compute empirical covariance
  Sigma_emp <- cov(estimates[, coords])

  # Create data frame
  df <- data.frame(x = x_est, y = y_est)

  # Plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
    geom_point(aes(x = x_true, y = y_true), color = "red", size = 4, shape = 3) +
    stat_ellipse(level = 0.95, color = "darkblue", linetype = "dashed") +
    labs(
      title = sprintf("Asymptotic Normality: Vertex %d", vertex_id),
      subtitle = sprintf("Coordinates %d vs %d (%d replicates)",
                         coords[1], coords[2], nrow(df)),
      x = sprintf("Coordinate %d", coords[1]),
      y = sprintf("Coordinate %d", coords[2])
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p)

  # Return empirical vs theoretical comparison
  list(
    plot = p,
    empirical_cov = Sigma_emp,
    theoretical_cov = Sigma_2d,
    empirical_mean = colMeans(estimates[, coords]),
    true_value = c(x_true, y_true)
  )
}

plot_marginal_histograms <- function(collected_data, vertex_idx = 1) {
  vertices <- collected_data$vertices_tracked
  vertex_id <- vertices[vertex_idx]
  X0 <- collected_data$X0_ref
  d <- dim(collected_data$estimates)[3]

  estimates <- collected_data$estimates[, vertex_idx, ]
  true_vals <- X0[vertex_id, ]

  # Theoretical SEs
  Sigma <- collected_data$theoretical_cov[vertex_idx, , ]
  theoretical_se <- sqrt(diag(Sigma))

  # Melt data for faceting
  df_list <- list()
  for (k in 1:d) {
    df_list[[k]] <- data.frame(
      coord = paste("Coordinate", k),
      value = estimates[, k],
      true_val = true_vals[k],
      theoretical_se = theoretical_se[k]
    )
  }
  df <- do.call(rbind, df_list)

  # Plot
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   fill = "steelblue", alpha = 0.6) +
    geom_vline(aes(xintercept = true_val), color = "red",
               linetype = "dashed", size = 1) +
    stat_function(fun = dnorm,
                  args = list(mean = df$true_val[1], sd = df$theoretical_se[1]),
                  color = "darkblue", size = 1) +
    facet_wrap(~ coord, scales = "free", ncol = 1) +
    labs(
      title = sprintf("Marginal Distributions: Vertex %d", vertex_id),
      subtitle = sprintf("%d replicates", nrow(estimates)),
      x = "Estimated Value",
      y = "Density"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p)
  return(p)
}

plot_qq_plots <- function(collected_data, vertex_idx = 1) {
  vertices <- collected_data$vertices_tracked
  vertex_id <- vertices[vertex_idx]
  X0 <- collected_data$X0_ref
  d <- dim(collected_data$estimates)[3]

  estimates <- collected_data$estimates[, vertex_idx, ]
  true_vals <- X0[vertex_id, ]

  # Theoretical SEs
  Sigma <- collected_data$theoretical_cov[vertex_idx, , ]
  theoretical_se <- sqrt(diag(Sigma))

  # Standardize
  df_list <- list()
  for (k in 1:d) {
    standardized <- (estimates[, k] - true_vals[k]) / theoretical_se[k]
    df_list[[k]] <- data.frame(
      coord = paste("Coordinate", k),
      standardized = standardized
    )
  }
  df <- do.call(rbind, df_list)

  # Q-Q plot
  p <- ggplot(df, aes(sample = standardized)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    facet_wrap(~ coord, ncol = 1) +
    labs(
      title = sprintf("Q-Q Plots: Vertex %d", vertex_id),
      subtitle = "Standardized estimates vs. N(0,1)",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p)
  return(p)
}

# ------------------------------------------------------------------------------
# RUN VISUALIZATION
# ------------------------------------------------------------------------------
cat("Collecting estimates across replicates...\n\n")

# Collect data
collected <- collect_estimates(
  n = 60,
  p_cov = 5,
  d = 3,
  n_reps = 500,
  vertices_to_track = c(10, 20, 30, 40, 50),
  seed = 456,
  verbose = TRUE
)

cat("\nCreating visualizations...\n")

# Create output directory for plots
dir.create("simulations/plots", showWarnings = FALSE, recursive = TRUE)

# 1. Scatter plots for different vertex/coordinate pairs
pdf("simulations/plots/scatter_asymptotic_normality.pdf", width = 10, height = 8)
for (v_idx in 1:length(collected$vertices_tracked)) {
  # Plot coordinates 1 vs 2
  result_12 <- plot_scatter_2d(collected, vertex_idx = v_idx, coords = c(1, 2))

  # Plot coordinates 1 vs 3
  result_13 <- plot_scatter_2d(collected, vertex_idx = v_idx, coords = c(1, 3))

  # Plot coordinates 2 vs 3
  result_23 <- plot_scatter_2d(collected, vertex_idx = v_idx, coords = c(2, 3))

  cat(sprintf("\nVertex %d:\n", collected$vertices_tracked[v_idx]))
  cat("Empirical cov (coords 1-2):\n")
  print(result_12$empirical_cov)
  cat("Theoretical cov (coords 1-2):\n")
  print(result_12$theoretical_cov)
}
dev.off()

# 2. Marginal histograms
pdf("simulations/plots/histograms_asymptotic_normality.pdf", width = 8, height = 10)
for (v_idx in 1:length(collected$vertices_tracked)) {
  plot_marginal_histograms(collected, vertex_idx = v_idx)
}
dev.off()

# 3. Q-Q plots
pdf("simulations/plots/qq_asymptotic_normality.pdf", width = 8, height = 10)
for (v_idx in 1:length(collected$vertices_tracked)) {
  plot_qq_plots(collected, vertex_idx = v_idx)
}
dev.off()

cat("\nPlots saved to simulations/plots/\n")
cat("  - scatter_asymptotic_normality.pdf\n")
cat("  - histograms_asymptotic_normality.pdf\n")
cat("  - qq_asymptotic_normality.pdf\n")

# Save collected data
saveRDS(collected, "simulations/asymptotic_normality_data.rds")
cat("\nData saved to simulations/asymptotic_normality_data.rds\n")
