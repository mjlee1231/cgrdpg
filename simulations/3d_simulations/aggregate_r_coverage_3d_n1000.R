#!/usr/bin/env Rscript
# Aggregate R coverage results from 100 replications (n=1000, p_cov=500)

cat("Aggregating 3D cgrdpg Coverage Results (R package, n=1000)\n")
cat("==========================================================\n\n")

results_dir <- "results_r_3d_coverage_n1000"
n_reps <- 100
n_vertices <- 1000

# Initialize storage
coverage_cgrdpg_true <- matrix(NA_real_, n_vertices, n_reps)
coverage_cgrdpg_plugin <- matrix(NA_real_, n_vertices, n_reps)
coverage_ase_true <- matrix(NA_real_, n_vertices, n_reps)
coverage_ase_plugin <- matrix(NA_real_, n_vertices, n_reps)

all_sse_cgrdpg <- numeric(n_reps)
all_sse_ase <- numeric(n_reps)
all_times_cgrdpg <- numeric(n_reps)
all_times_ase <- numeric(n_reps)
all_converged <- logical(n_reps)
all_iters <- numeric(n_reps)

# Load each replication
n_loaded <- 0
for (rep in 1:n_reps) {
  filename <- file.path(results_dir, sprintf("rep_%03d.rds", rep))

  if (file.exists(filename)) {
    data <- readRDS(filename)

    # Store vertex-wise coverage
    coverage_cgrdpg_true[, rep] <- data$results_mat[, "cgrdpg_true"]
    coverage_cgrdpg_plugin[, rep] <- data$results_mat[, "cgrdpg_plugin"]
    coverage_ase_true[, rep] <- data$results_mat[, "ase_true"]
    coverage_ase_plugin[, rep] <- data$results_mat[, "ase_plugin"]

    # Store summary statistics
    all_sse_cgrdpg[rep] <- data$sse_cgrdpg
    all_sse_ase[rep] <- data$sse_ase
    all_times_cgrdpg[rep] <- data$cgrdpg_time
    all_times_ase[rep] <- data$ase_time
    all_converged[rep] <- data$converged
    all_iters[rep] <- data$iters

    n_loaded <- n_loaded + 1
  } else {
    cat(sprintf("Warning: Rep %d not found\n", rep))
  }
}

cat(sprintf("Loaded %d/%d replications\n\n", n_loaded, n_reps))

# Compute vertex-wise coverage rates (across replications)
vertexwise_cov_cgrdpg_true <- rowMeans(coverage_cgrdpg_true, na.rm = TRUE)
vertexwise_cov_cgrdpg_plugin <- rowMeans(coverage_cgrdpg_plugin, na.rm = TRUE)
vertexwise_cov_ase_true <- rowMeans(coverage_ase_true, na.rm = TRUE)
vertexwise_cov_ase_plugin <- rowMeans(coverage_ase_plugin, na.rm = TRUE)

# Overall coverage (average across all vertices)
overall_cov_cgrdpg_true <- mean(vertexwise_cov_cgrdpg_true, na.rm = TRUE)
overall_cov_cgrdpg_plugin <- mean(vertexwise_cov_cgrdpg_plugin, na.rm = TRUE)
overall_cov_ase_true <- mean(vertexwise_cov_ase_true, na.rm = TRUE)
overall_cov_ase_plugin <- mean(vertexwise_cov_ase_plugin, na.rm = TRUE)

# Summary
cat("========================================\n")
cat("Vertex-wise Coverage Summary (n=1000)\n")
cat("========================================\n\n")

cat("cgrdpg-TRUE:\n")
cat(sprintf("  Overall (avg across vertices): %.2f%%\n", 100 * overall_cov_cgrdpg_true))
cat(sprintf("  Vertex coverage range:         [%.2f%%, %.2f%%]\n",
            100 * min(vertexwise_cov_cgrdpg_true, na.rm = TRUE),
            100 * max(vertexwise_cov_cgrdpg_true, na.rm = TRUE)))
cat(sprintf("  Std across vertices:           %.2f%%\n\n",
            100 * sd(vertexwise_cov_cgrdpg_true, na.rm = TRUE)))

cat("cgrdpg-PLUGIN:\n")
cat(sprintf("  Overall (avg across vertices): %.2f%%\n", 100 * overall_cov_cgrdpg_plugin))
cat(sprintf("  Vertex coverage range:         [%.2f%%, %.2f%%]\n",
            100 * min(vertexwise_cov_cgrdpg_plugin, na.rm = TRUE),
            100 * max(vertexwise_cov_cgrdpg_plugin, na.rm = TRUE)))
cat(sprintf("  Std across vertices:           %.2f%%\n\n",
            100 * sd(vertexwise_cov_cgrdpg_plugin, na.rm = TRUE)))

cat("ASE-TRUE:\n")
cat(sprintf("  Overall (avg across vertices): %.2f%%\n", 100 * overall_cov_ase_true))
cat(sprintf("  Vertex coverage range:         [%.2f%%, %.2f%%]\n",
            100 * min(vertexwise_cov_ase_true, na.rm = TRUE),
            100 * max(vertexwise_cov_ase_true, na.rm = TRUE)))
cat(sprintf("  Std across vertices:           %.2f%%\n\n",
            100 * sd(vertexwise_cov_ase_true, na.rm = TRUE)))

cat("ASE-PLUGIN:\n")
cat(sprintf("  Overall (avg across vertices): %.2f%%\n", 100 * overall_cov_ase_plugin))
cat(sprintf("  Vertex coverage range:         [%.2f%%, %.2f%%]\n",
            100 * min(vertexwise_cov_ase_plugin, na.rm = TRUE),
            100 * max(vertexwise_cov_ase_plugin, na.rm = TRUE)))
cat(sprintf("  Std across vertices:           %.2f%%\n\n",
            100 * sd(vertexwise_cov_ase_plugin, na.rm = TRUE)))

cat("Optimization:\n")
cat(sprintf("  cgrdpg Convergence: %.1f%%\n", 100 * mean(all_converged, na.rm = TRUE)))
cat(sprintf("  cgrdpg Mean iters:  %.1f\n", mean(all_iters, na.rm = TRUE)))
cat(sprintf("  cgrdpg Mean time:   %.1f sec\n", mean(all_times_cgrdpg, na.rm = TRUE)))
cat(sprintf("  ASE Mean time:      %.1f sec\n\n", mean(all_times_ase, na.rm = TRUE)))

cat("SSE Distribution (cgrdpg):\n")
cat(sprintf("  Mean:   %.4f\n", mean(all_sse_cgrdpg, na.rm = TRUE)))
cat(sprintf("  Std:    %.4f\n", sd(all_sse_cgrdpg, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(all_sse_cgrdpg, na.rm = TRUE)))
cat(sprintf("  Range:  [%.4f, %.4f]\n\n", min(all_sse_cgrdpg), max(all_sse_cgrdpg)))

cat("SSE Distribution (ASE):\n")
cat(sprintf("  Mean:   %.4f\n", mean(all_sse_ase, na.rm = TRUE)))
cat(sprintf("  Std:    %.4f\n", sd(all_sse_ase, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(all_sse_ase, na.rm = TRUE)))
cat(sprintf("  Range:  [%.4f, %.4f]\n\n", min(all_sse_ase), max(all_sse_ase)))

# Save aggregated results
results <- list(
  n_reps = n_reps,
  n_loaded = n_loaded,
  n_vertices = n_vertices,
  vertexwise_cov_cgrdpg_true = vertexwise_cov_cgrdpg_true,
  vertexwise_cov_cgrdpg_plugin = vertexwise_cov_cgrdpg_plugin,
  vertexwise_cov_ase_true = vertexwise_cov_ase_true,
  vertexwise_cov_ase_plugin = vertexwise_cov_ase_plugin,
  overall_cov_cgrdpg_true = overall_cov_cgrdpg_true,
  overall_cov_cgrdpg_plugin = overall_cov_cgrdpg_plugin,
  overall_cov_ase_true = overall_cov_ase_true,
  overall_cov_ase_plugin = overall_cov_ase_plugin,
  all_sse_cgrdpg = all_sse_cgrdpg,
  all_sse_ase = all_sse_ase,
  all_times_cgrdpg = all_times_cgrdpg,
  all_times_ase = all_times_ase,
  all_converged = all_converged,
  all_iters = all_iters
)

saveRDS(results, file.path(results_dir, "aggregated_results.rds"))

cat(sprintf("\nAggregated results saved to: %s\n",
            file.path(results_dir, "aggregated_results.rds")))
