#!/usr/bin/env Rscript
# Aggregate comprehensive 1D simulation results across 100 replications
# Computes vertex-wise coverage for Fisher vs ASE vs OSE

cat("============================================================================\n")
cat("  AGGREGATING 1D COMPREHENSIVE SIMULATION RESULTS\n")
cat("  Computing vertex-wise coverage across 100 replications\n")
cat("============================================================================\n\n")

results_dir <- "comprehensive_1d_n500_results"

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop(sprintf("Results directory not found: %s\n", results_dir))
}

# Find all result files
result_files <- list.files(results_dir, pattern = "rep_\\d{3}\\.rds", full.names = TRUE)
n_reps <- length(result_files)

if (n_reps == 0) {
  stop(sprintf("No result files found in %s\n", results_dir))
}

cat(sprintf("Found %d replications\n", n_reps))

# Read first file to get dimensions
first_result <- readRDS(result_files[1])
n <- first_result$n
cat(sprintf("Number of vertices: %d\n\n", n))

# Initialize storage
fisher_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
fisher_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ase_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ase_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)

# Storage for SSE and timing
sse_fisher <- numeric(n_reps)
sse_ase <- numeric(n_reps)
sse_ose <- numeric(n_reps)

time_fisher <- numeric(n_reps)
time_ase <- numeric(n_reps)
time_ose <- numeric(n_reps)
time_coverage <- numeric(n_reps)

cat("Reading all replications...\n")
for (r in 1:n_reps) {
  if (r %% 10 == 0) cat(sprintf("  Processing replication %d/%d\n", r, n_reps))

  result <- readRDS(result_files[r])

  # Coverage
  fisher_true_matrix[, r] <- result$coverage$fisher_true
  fisher_plugin_matrix[, r] <- result$coverage$fisher_plugin
  ase_true_matrix[, r] <- result$coverage$ase_true
  ase_plugin_matrix[, r] <- result$coverage$ase_plugin
  ose_true_matrix[, r] <- result$coverage$ose_true
  ose_plugin_matrix[, r] <- result$coverage$ose_plugin

  # SSE
  sse_fisher[r] <- result$sse$fisher
  sse_ase[r] <- result$sse$ase
  sse_ose[r] <- result$sse$ose

  # Timing
  time_fisher[r] <- result$timing$fisher
  time_ase[r] <- result$timing$ase
  time_ose[r] <- result$timing$ose
  time_coverage[r] <- result$timing$coverage
}

cat("\nComputing vertex-wise coverage...\n")

# Vertex-wise coverage (proportion of reps where vertex was covered)
vertex_fisher_true <- rowMeans(fisher_true_matrix, na.rm = TRUE)
vertex_fisher_plugin <- rowMeans(fisher_plugin_matrix, na.rm = TRUE)
vertex_ase_true <- rowMeans(ase_true_matrix, na.rm = TRUE)
vertex_ase_plugin <- rowMeans(ase_plugin_matrix, na.rm = TRUE)
vertex_ose_true <- rowMeans(ose_true_matrix, na.rm = TRUE)
vertex_ose_plugin <- rowMeans(ose_plugin_matrix, na.rm = TRUE)

# Overall coverage (average across all vertices and reps)
overall_fisher_true <- mean(fisher_true_matrix, na.rm = TRUE)
overall_fisher_plugin <- mean(fisher_plugin_matrix, na.rm = TRUE)
overall_ase_true <- mean(ase_true_matrix, na.rm = TRUE)
overall_ase_plugin <- mean(ase_plugin_matrix, na.rm = TRUE)
overall_ose_true <- mean(ose_true_matrix, na.rm = TRUE)
overall_ose_plugin <- mean(ose_plugin_matrix, na.rm = TRUE)

# Print results
cat("\n============================================================================\n")
cat("  OVERALL COVERAGE RESULTS (n=500, 100 replications)\n")
cat("============================================================================\n\n")

cat("FISHER-SCORING (with covariates):\n")
cat(sprintf("  TRUE precision:   %.2f%%\n", 100 * overall_fisher_true))
cat(sprintf("  PLUGIN precision: %.2f%%\n\n", 100 * overall_fisher_plugin))

cat("ASE (no covariates):\n")
cat(sprintf("  TRUE precision:   %.2f%%\n", 100 * overall_ase_true))
cat(sprintf("  PLUGIN precision: %.2f%%\n\n", 100 * overall_ase_plugin))

cat("OSE (no covariates):\n")
cat(sprintf("  TRUE precision:   %.2f%%\n", 100 * overall_ose_true))
cat(sprintf("  PLUGIN precision: %.2f%%\n\n", 100 * overall_ose_plugin))

cat("============================================================================\n")
cat("  SSE COMPARISON\n")
cat("============================================================================\n\n")

cat(sprintf("FISHER: Mean=%.4f, SD=%.4f, Min=%.4f, Max=%.4f\n",
            mean(sse_fisher), sd(sse_fisher), min(sse_fisher), max(sse_fisher)))
cat(sprintf("ASE:    Mean=%.4f, SD=%.4f, Min=%.4f, Max=%.4f\n",
            mean(sse_ase), sd(sse_ase), min(sse_ase), max(sse_ase)))
cat(sprintf("OSE:    Mean=%.4f, SD=%.4f, Min=%.4f, Max=%.4f\n\n",
            mean(sse_ose), sd(sse_ose), min(sse_ose), max(sse_ose)))

cat("============================================================================\n")
cat("  TIMING COMPARISON\n")
cat("============================================================================\n\n")

cat(sprintf("FISHER:   Mean=%.2fs, SD=%.2fs\n", mean(time_fisher), sd(time_fisher)))
cat(sprintf("ASE:      Mean=%.2fs, SD=%.2fs\n", mean(time_ase), sd(time_ase)))
cat(sprintf("OSE:      Mean=%.2fs, SD=%.2fs\n", mean(time_ose), sd(time_ose)))
cat(sprintf("Coverage: Mean=%.2fs, SD=%.2fs\n", mean(time_coverage), sd(time_coverage)))
cat(sprintf("Total:    Mean=%.2f min per replication\n\n",
            mean(time_fisher + time_ase + time_ose + time_coverage)/60))

# Save aggregated results
aggregated <- list(
  n = n,
  n_reps = n_reps,
  overall_coverage = list(
    fisher_true = overall_fisher_true,
    fisher_plugin = overall_fisher_plugin,
    ase_true = overall_ase_true,
    ase_plugin = overall_ase_plugin,
    ose_true = overall_ose_true,
    ose_plugin = overall_ose_plugin
  ),
  vertex_coverage = list(
    fisher_true = vertex_fisher_true,
    fisher_plugin = vertex_fisher_plugin,
    ase_true = vertex_ase_true,
    ase_plugin = vertex_ase_plugin,
    ose_true = vertex_ose_true,
    ose_plugin = vertex_ose_plugin
  ),
  coverage_matrices = list(
    fisher_true = fisher_true_matrix,
    fisher_plugin = fisher_plugin_matrix,
    ase_true = ase_true_matrix,
    ase_plugin = ase_plugin_matrix,
    ose_true = ose_true_matrix,
    ose_plugin = ose_plugin_matrix
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
    coverage = time_coverage
  )
)

saveRDS(aggregated, "comprehensive_1d_n500_aggregated.rds")
cat("Aggregated results saved to: comprehensive_1d_n500_aggregated.rds\n")

# Save CSV
df <- data.frame(
  vertex = 1:n,
  fisher_true = vertex_fisher_true,
  fisher_plugin = vertex_fisher_plugin,
  ase_true = vertex_ase_true,
  ase_plugin = vertex_ase_plugin,
  ose_true = vertex_ose_true,
  ose_plugin = vertex_ose_plugin
)
write.csv(df, "comprehensive_1d_n500_vertex_coverage.csv", row.names = FALSE)
cat("Vertex-wise CSV saved to: comprehensive_1d_n500_vertex_coverage.csv\n")

cat("\n============================================================================\n")
cat("  AGGREGATION COMPLETE\n")
cat("============================================================================\n")
