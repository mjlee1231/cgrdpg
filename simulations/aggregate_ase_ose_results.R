#!/usr/bin/env Rscript
# Aggregate ASE vs OSE results across 100 replications
# Computes vertex-wise coverage for both methods with TRUE and PLUGIN G_in

cat("============================================================================\n")
cat("  AGGREGATING ASE vs OSE RESULTS\n")
cat("  100 Replications, n=500\n")
cat("============================================================================\n\n")

# Find all result files
result_dir <- "ase_ose_results_n500"
if (!dir.exists(result_dir)) {
  stop(sprintf("Directory '%s' not found! Make sure simulations have completed.", result_dir))
}

result_files <- list.files(path = result_dir,
                           pattern = "^ase_ose_rep[0-9]{3}\\.rds$",
                           full.names = TRUE)

if (length(result_files) == 0) {
  stop(sprintf("No result files found in %s!", result_dir))
}

cat(sprintf("Found %d result files\n", length(result_files)))

if (length(result_files) < 100) {
  cat(sprintf("WARNING: Expected 100 files, found only %d\n", length(result_files)))
  cat("Proceeding with available files...\n\n")
} else {
  cat("All 100 replications found!\n\n")
}

# Read all results
cat("Reading result files...\n")
all_results <- lapply(result_files, readRDS)

# Extract parameters
n <- all_results[[1]]$n
d <- all_results[[1]]$d
n_reps <- length(all_results)

cat(sprintf("Parameters: n=%d, d=%d, n_reps=%d\n\n", n, d, n_reps))

# Initialize matrices for per-vertex coverage across reps (n × n_reps)
ase_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ase_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
ose_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)

# Extract per-replication statistics
rep_ids <- integer(n_reps)
rep_sse_ase <- numeric(n_reps)
rep_sse_ose <- numeric(n_reps)
rep_cov_ase_true <- numeric(n_reps)
rep_cov_ase_plugin <- numeric(n_reps)
rep_cov_ose_true <- numeric(n_reps)
rep_cov_ose_plugin <- numeric(n_reps)

cat("Extracting per-vertex indicators from each replication...\n")
for (j in 1:n_reps) {
  rep_ids[j] <- all_results[[j]]$rep_id
  rep_sse_ase[j] <- all_results[[j]]$sse$ase
  rep_sse_ose[j] <- all_results[[j]]$sse$ose
  rep_cov_ase_true[j] <- all_results[[j]]$overall_coverage$ase_true
  rep_cov_ase_plugin[j] <- all_results[[j]]$overall_coverage$ase_plugin
  rep_cov_ose_true[j] <- all_results[[j]]$overall_coverage$ose_true
  rep_cov_ose_plugin[j] <- all_results[[j]]$overall_coverage$ose_plugin

  # Store per-vertex coverage indicators
  ase_true_matrix[, j] <- all_results[[j]]$coverage_df$ase_true
  ase_plugin_matrix[, j] <- all_results[[j]]$coverage_df$ase_plugin
  ose_true_matrix[, j] <- all_results[[j]]$coverage_df$ose_true
  ose_plugin_matrix[, j] <- all_results[[j]]$coverage_df$ose_plugin
}

cat("Computing vertex-wise coverage rates...\n\n")

# Compute vertex-wise coverage (average across replications for each vertex)
vertex_cov_ase_true <- apply(ase_true_matrix, 1, mean)
vertex_cov_ase_plugin <- apply(ase_plugin_matrix, 1, mean)
vertex_cov_ose_true <- apply(ose_true_matrix, 1, mean)
vertex_cov_ose_plugin <- apply(ose_plugin_matrix, 1, mean)

cat("============================================================================\n")
cat("  VERTEX-WISE COVERAGE RESULTS\n")
cat("============================================================================\n\n")

cat("--- ASE with TRUE G_in ---\n")
cat(sprintf("  Mean:   %.2f%%\n", 100*mean(vertex_cov_ase_true)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_cov_ase_true)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_cov_ase_true)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_cov_ase_true), 100*max(vertex_cov_ase_true)))

cat("--- ASE with PLUGIN G_in ---\n")
cat(sprintf("  Mean:   %.2f%%\n", 100*mean(vertex_cov_ase_plugin)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_cov_ase_plugin)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_cov_ase_plugin)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_cov_ase_plugin), 100*max(vertex_cov_ase_plugin)))

cat("--- OSE with TRUE G_in ---\n")
cat(sprintf("  Mean:   %.2f%%\n", 100*mean(vertex_cov_ose_true)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_cov_ose_true)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_cov_ose_true)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_cov_ose_true), 100*max(vertex_cov_ose_true)))

cat("--- OSE with PLUGIN G_in ---\n")
cat(sprintf("  Mean:   %.2f%%\n", 100*mean(vertex_cov_ose_plugin)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_cov_ose_plugin)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_cov_ose_plugin)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_cov_ose_plugin), 100*max(vertex_cov_ose_plugin)))

# Comparisons
cat("--- Method Comparisons ---\n")
cat(sprintf("OSE vs ASE (TRUE G_in):   %.2f%% improvement\n",
            100*(mean(vertex_cov_ose_true) - mean(vertex_cov_ase_true))))
cat(sprintf("OSE vs ASE (PLUGIN G_in): %.2f%% improvement\n",
            100*(mean(vertex_cov_ose_plugin) - mean(vertex_cov_ase_plugin))))
cat(sprintf("PLUGIN vs TRUE (ASE):     %.2f%% difference\n",
            100*(mean(vertex_cov_ase_plugin) - mean(vertex_cov_ase_true))))
cat(sprintf("PLUGIN vs TRUE (OSE):     %.2f%% difference\n\n",
            100*(mean(vertex_cov_ose_plugin) - mean(vertex_cov_ose_true))))

# Example vertices
cat("--- Example Vertices ---\n")
example_vertices <- c(1, 100, 250, 400, 500)
cat(sprintf("%-8s %-12s %-12s %-12s %-12s\n",
            "Vertex", "ASE-TRUE%", "ASE-PLUG%", "OSE-TRUE%", "OSE-PLUG%"))
cat(paste(rep("-", 65), collapse=""), "\n")
for (v in example_vertices) {
  cat(sprintf("%-8d %-12.1f %-12.1f %-12.1f %-12.1f\n",
              v,
              100*vertex_cov_ase_true[v],
              100*vertex_cov_ase_plugin[v],
              100*vertex_cov_ose_true[v],
              100*vertex_cov_ose_plugin[v]))
}

cat("\n============================================================================\n")
cat("  PER-REPLICATION SUMMARY\n")
cat("============================================================================\n\n")

cat(sprintf("SSE (ASE):  Mean=%.4f, SD=%.4f\n", mean(rep_sse_ase), sd(rep_sse_ase)))
cat(sprintf("SSE (OSE):  Mean=%.4f, SD=%.4f\n\n", mean(rep_sse_ose), sd(rep_sse_ose)))

cat(sprintf("Coverage ASE-TRUE:   Mean=%.1f%%, SD=%.1f%%\n",
            100*mean(rep_cov_ase_true), 100*sd(rep_cov_ase_true)))
cat(sprintf("Coverage ASE-PLUGIN: Mean=%.1f%%, SD=%.1f%%\n",
            100*mean(rep_cov_ase_plugin), 100*sd(rep_cov_ase_plugin)))
cat(sprintf("Coverage OSE-TRUE:   Mean=%.1f%%, SD=%.1f%%\n",
            100*mean(rep_cov_ose_true), 100*sd(rep_cov_ose_true)))
cat(sprintf("Coverage OSE-PLUGIN: Mean=%.1f%%, SD=%.1f%%\n\n",
            100*mean(rep_cov_ose_plugin), 100*sd(rep_cov_ose_plugin)))

# Save aggregated results
aggregated_results <- list(
  n = n,
  d = d,
  n_reps = n_reps,

  # Vertex-wise coverage (length n = 500)
  vertex_cov_ase_true = vertex_cov_ase_true,
  vertex_cov_ase_plugin = vertex_cov_ase_plugin,
  vertex_cov_ose_true = vertex_cov_ose_true,
  vertex_cov_ose_plugin = vertex_cov_ose_plugin,

  # Per-vertex matrices (n × n_reps)
  ase_true_matrix = ase_true_matrix,
  ase_plugin_matrix = ase_plugin_matrix,
  ose_true_matrix = ose_true_matrix,
  ose_plugin_matrix = ose_plugin_matrix,

  # Per-replication statistics
  rep_ids = rep_ids,
  rep_sse_ase = rep_sse_ase,
  rep_sse_ose = rep_sse_ose,
  rep_cov_ase_true = rep_cov_ase_true,
  rep_cov_ase_plugin = rep_cov_ase_plugin,
  rep_cov_ose_true = rep_cov_ose_true,
  rep_cov_ose_plugin = rep_cov_ose_plugin
)

output_file <- file.path(result_dir, "ase_ose_aggregated.rds")
saveRDS(aggregated_results, output_file)
cat(sprintf("Saved: %s\n\n", output_file))

# Also save as CSV
csv_data <- data.frame(
  vertex = 1:n,
  ase_true = vertex_cov_ase_true,
  ase_plugin = vertex_cov_ase_plugin,
  ose_true = vertex_cov_ose_true,
  ose_plugin = vertex_cov_ose_plugin
)
csv_file <- file.path(result_dir, "ase_ose_vertex_wise_coverage.csv")
write.csv(csv_data, csv_file, row.names = FALSE)
cat(sprintf("Saved: %s\n\n", csv_file))

cat("============================================================================\n")
cat("  AGGREGATION COMPLETE!\n")
cat("============================================================================\n")
