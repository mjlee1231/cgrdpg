#!/usr/bin/env Rscript
# Compute vertex-wise coverage across replications for ASE vs OSE (n=3000)
# For each vertex i: what is its coverage rate across 100 replications?

cat("============================================================================\n")
cat("  VERTEX-WISE COVERAGE ACROSS REPLICATIONS: ASE vs OSE (n=3000)\n")
cat("  Computing coverage for each vertex across 100 replications\n")
cat("============================================================================\n\n")

results_dir <- "ase_ose_results_n3000"

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop(sprintf("Results directory not found: %s\n", results_dir))
}

# Find all result files
result_files <- list.files(results_dir, pattern = "ase_ose_n3000_rep\\d{3}\\.rds", full.names = TRUE)
n_reps <- length(result_files)

if (n_reps == 0) {
  stop(sprintf("No result files found in %s\n", results_dir))
}

cat(sprintf("Found %d replications\n", n_reps))

# Read first file to get n
first_result <- readRDS(result_files[1])
n <- first_result$n
cat(sprintf("Number of vertices: %d\n\n", n))

# Initialize matrices to store coverage indicators
# Rows = vertices, Columns = replications
ase_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ase_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)

# Also collect SSE and timing
sse_ase_vec <- numeric(n_reps)
sse_ose_vec <- numeric(n_reps)
time_ase_vec <- numeric(n_reps)
time_ose_vec <- numeric(n_reps)
time_coverage_vec <- numeric(n_reps)

cat("Reading all replications...\n")
for (r in 1:n_reps) {
  if (r %% 10 == 0) cat(sprintf("  Processing replication %d/%d\n", r, n_reps))

  result <- readRDS(result_files[r])

  # Store coverage indicators (TRUE/FALSE for each vertex)
  ase_true_matrix[, r] <- result$coverage_df$ase_true
  ase_plugin_matrix[, r] <- result$coverage_df$ase_plugin
  ose_true_matrix[, r] <- result$coverage_df$ose_true
  ose_plugin_matrix[, r] <- result$coverage_df$ose_plugin

  # Store SSE and timing
  sse_ase_vec[r] <- result$sse$ase
  sse_ose_vec[r] <- result$sse$ose
  time_ase_vec[r] <- result$timing$ase_time
  time_ose_vec[r] <- result$timing$ose_time
  time_coverage_vec[r] <- result$timing$coverage_time
}

cat("\nComputing vertex-wise coverage...\n")

# For each vertex, compute coverage across replications
# (proportion of replications where vertex was covered)
vertex_coverage_ase_true <- rowMeans(ase_true_matrix, na.rm = TRUE)
vertex_coverage_ase_plugin <- rowMeans(ase_plugin_matrix, na.rm = TRUE)
vertex_coverage_ose_true <- rowMeans(ose_true_matrix, na.rm = TRUE)
vertex_coverage_ose_plugin <- rowMeans(ose_plugin_matrix, na.rm = TRUE)

# Summary statistics
cat("\n============================================================================\n")
cat("  VERTEX-WISE COVERAGE SUMMARY (n=3000)\n")
cat("============================================================================\n\n")

cat("ASE-TRUE coverage (per vertex across replications):\n")
cat(sprintf("  Mean:   %.4f (%.2f%%)\n", mean(vertex_coverage_ase_true, na.rm=TRUE),
            100*mean(vertex_coverage_ase_true, na.rm=TRUE)))
cat(sprintf("  Median: %.4f\n", median(vertex_coverage_ase_true, na.rm=TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(vertex_coverage_ase_true, na.rm=TRUE)))
cat(sprintf("  Min:    %.4f (%.2f%%)\n", min(vertex_coverage_ase_true, na.rm=TRUE),
            100*min(vertex_coverage_ase_true, na.rm=TRUE)))
cat(sprintf("  Max:    %.4f (%.2f%%)\n\n", max(vertex_coverage_ase_true, na.rm=TRUE),
            100*max(vertex_coverage_ase_true, na.rm=TRUE)))

cat("ASE-PLUGIN coverage (per vertex across replications):\n")
cat(sprintf("  Mean:   %.4f (%.2f%%)\n", mean(vertex_coverage_ase_plugin, na.rm=TRUE),
            100*mean(vertex_coverage_ase_plugin, na.rm=TRUE)))
cat(sprintf("  Median: %.4f\n", median(vertex_coverage_ase_plugin, na.rm=TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(vertex_coverage_ase_plugin, na.rm=TRUE)))
cat(sprintf("  Min:    %.4f (%.2f%%)\n", min(vertex_coverage_ase_plugin, na.rm=TRUE),
            100*min(vertex_coverage_ase_plugin, na.rm=TRUE)))
cat(sprintf("  Max:    %.4f (%.2f%%)\n\n", max(vertex_coverage_ase_plugin, na.rm=TRUE),
            100*max(vertex_coverage_ase_plugin, na.rm=TRUE)))

cat("OSE-TRUE coverage (per vertex across replications):\n")
cat(sprintf("  Mean:   %.4f (%.2f%%)\n", mean(vertex_coverage_ose_true, na.rm=TRUE),
            100*mean(vertex_coverage_ose_true, na.rm=TRUE)))
cat(sprintf("  Median: %.4f\n", median(vertex_coverage_ose_true, na.rm=TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(vertex_coverage_ose_true, na.rm=TRUE)))
cat(sprintf("  Min:    %.4f (%.2f%%)\n", min(vertex_coverage_ose_true, na.rm=TRUE),
            100*min(vertex_coverage_ose_true, na.rm=TRUE)))
cat(sprintf("  Max:    %.4f (%.2f%%)\n\n", max(vertex_coverage_ose_true, na.rm=TRUE),
            100*max(vertex_coverage_ose_true, na.rm=TRUE)))

cat("OSE-PLUGIN coverage (per vertex across replications):\n")
cat(sprintf("  Mean:   %.4f (%.2f%%)\n", mean(vertex_coverage_ose_plugin, na.rm=TRUE),
            100*mean(vertex_coverage_ose_plugin, na.rm=TRUE)))
cat(sprintf("  Median: %.4f\n", median(vertex_coverage_ose_plugin, na.rm=TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(vertex_coverage_ose_plugin, na.rm=TRUE)))
cat(sprintf("  Min:    %.4f (%.2f%%)\n", min(vertex_coverage_ose_plugin, na.rm=TRUE),
            100*min(vertex_coverage_ose_plugin, na.rm=TRUE)))
cat(sprintf("  Max:    %.4f (%.2f%%)\n\n", max(vertex_coverage_ose_plugin, na.rm=TRUE),
            100*max(vertex_coverage_ose_plugin, na.rm=TRUE)))

# SSE summary
cat("============================================================================\n")
cat("  SSE SUMMARY\n")
cat("============================================================================\n\n")

cat("ASE SSE:\n")
cat(sprintf("  Mean: %.4f, SD: %.4f\n", mean(sse_ase_vec), sd(sse_ase_vec)))
cat(sprintf("  Min:  %.4f, Max: %.4f\n\n", min(sse_ase_vec), max(sse_ase_vec)))

cat("OSE SSE:\n")
cat(sprintf("  Mean: %.4f, SD: %.4f\n", mean(sse_ose_vec), sd(sse_ose_vec)))
cat(sprintf("  Min:  %.4f, Max: %.4f\n\n", min(sse_ose_vec), max(sse_ose_vec)))

cat(sprintf("SSE Improvement (OSE vs ASE): %.2f%%\n\n",
            100 * (mean(sse_ase_vec) - mean(sse_ose_vec)) / mean(sse_ase_vec)))

# Timing summary
cat("============================================================================\n")
cat("  TIMING SUMMARY\n")
cat("============================================================================\n\n")

cat(sprintf("ASE time:      Mean=%.1fs, SD=%.1fs\n", mean(time_ase_vec), sd(time_ase_vec)))
cat(sprintf("OSE time:      Mean=%.1fs, SD=%.1fs\n", mean(time_ose_vec), sd(time_ose_vec)))
cat(sprintf("Coverage time: Mean=%.1fs, SD=%.1fs\n", mean(time_coverage_vec), sd(time_coverage_vec)))
cat(sprintf("Total time:    Mean=%.1f min per replication\n\n",
            mean(time_ase_vec + time_ose_vec + time_coverage_vec)/60))

# Save results
vertexwise_results <- list(
  n = n,
  n_reps = n_reps,
  vertex_coverage = list(
    ase_true = vertex_coverage_ase_true,
    ase_plugin = vertex_coverage_ase_plugin,
    ose_true = vertex_coverage_ose_true,
    ose_plugin = vertex_coverage_ose_plugin
  ),
  coverage_matrices = list(
    ase_true = ase_true_matrix,
    ase_plugin = ase_plugin_matrix,
    ose_true = ose_true_matrix,
    ose_plugin = ose_plugin_matrix
  ),
  sse = list(
    ase = sse_ase_vec,
    ose = sse_ose_vec
  ),
  timing = list(
    ase = time_ase_vec,
    ose = time_ose_vec,
    coverage = time_coverage_vec
  )
)

saveRDS(vertexwise_results, "vertexwise_coverage_ase_ose_n3000.rds")
cat("Vertex-wise coverage saved to: vertexwise_coverage_ase_ose_n3000.rds\n")

# Save CSV
df <- data.frame(
  vertex = 1:n,
  ase_true = vertex_coverage_ase_true,
  ase_plugin = vertex_coverage_ase_plugin,
  ose_true = vertex_coverage_ose_true,
  ose_plugin = vertex_coverage_ose_plugin
)
write.csv(df, "vertexwise_coverage_ase_ose_n3000.csv", row.names = FALSE)
cat("Vertex-wise coverage CSV saved to: vertexwise_coverage_ase_ose_n3000.csv\n")

cat("\n============================================================================\n")
cat("  VERTEX-WISE AGGREGATION COMPLETE\n")
cat("============================================================================\n")
