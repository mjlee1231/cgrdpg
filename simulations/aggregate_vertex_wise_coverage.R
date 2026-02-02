#!/usr/bin/env Rscript
# Aggregate vertex-wise coverage across 100 replications
# Computes coverage rate for EACH vertex across all reps

cat("============================================================================\n")
cat("  AGGREGATING VERTEX-WISE COVERAGE RESULTS\n")
cat("  2D SPREAD - 100 Replications\n")
cat("============================================================================\n\n")

# Find all result files in the vertex_wise_results_n500_optimized directory
result_dir <- "vertex_wise_results_n500_optimized"
if (!dir.exists(result_dir)) {
  stop(sprintf("Directory '%s' not found! Make sure simulations have completed.", result_dir))
}

result_files <- list.files(path = result_dir,
                           pattern = "^vertex_wise_spread_rep[0-9]{3}\\.rds$",
                           full.names = TRUE)

if (length(result_files) == 0) {
  stop(sprintf("No result files found in %s! Expected: vertex_wise_spread_rep001.rds through rep100.rds", result_dir))
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

# Extract parameters from first result
n <- all_results[[1]]$n
p_cov <- all_results[[1]]$p_cov
d <- all_results[[1]]$d
n_reps <- length(all_results)

cat(sprintf("Parameters: n=%d vertices, p_cov=%d, d=%d, n_reps=%d\n\n", n, p_cov, d, n_reps))

# Initialize matrices to store per-vertex indicators across all reps
# Rows = vertices (1 to n), Columns = replications
in_ellipse_true_matrix <- matrix(NA, nrow = n, ncol = n_reps)
in_ellipse_plugin_matrix <- matrix(NA, nrow = n, ncol = n_reps)

# Extract per-replication statistics
rep_ids <- integer(n_reps)
rep_SSEs <- numeric(n_reps)
rep_fit_times <- numeric(n_reps)
rep_coverage_true <- numeric(n_reps)
rep_coverage_plugin <- numeric(n_reps)

cat("Extracting per-vertex indicators from each replication...\n")
for (j in 1:n_reps) {
  rep_ids[j] <- all_results[[j]]$rep
  rep_SSEs[j] <- all_results[[j]]$SSE
  rep_fit_times[j] <- all_results[[j]]$fit_time
  rep_coverage_true[j] <- all_results[[j]]$coverage_true
  rep_coverage_plugin[j] <- all_results[[j]]$coverage_plugin

  # Store per-vertex indicators
  in_ellipse_true_matrix[, j] <- all_results[[j]]$in_ellipse_true
  in_ellipse_plugin_matrix[, j] <- all_results[[j]]$in_ellipse_plugin
}

cat("Computing vertex-wise coverage rates...\n\n")

# Compute vertex-wise coverage rate (average across replications for each vertex)
vertex_coverage_true <- apply(in_ellipse_true_matrix, 1, function(x) mean(x, na.rm = TRUE))
vertex_coverage_plugin <- apply(in_ellipse_plugin_matrix, 1, function(x) mean(x, na.rm = TRUE))

# Count NAs per vertex
vertex_n_na_true <- apply(in_ellipse_true_matrix, 1, function(x) sum(is.na(x)))
vertex_n_na_plugin <- apply(in_ellipse_plugin_matrix, 1, function(x) sum(is.na(x)))

cat("============================================================================\n")
cat("  VERTEX-WISE COVERAGE RESULTS\n")
cat("============================================================================\n\n")

cat("--- Vertex-wise Coverage Statistics ---\n")
cat(sprintf("TRUE Coverage:\n"))
cat(sprintf("  Mean:   %.2f%% (average across all vertices)\n", 100*mean(vertex_coverage_true)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_coverage_true)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_coverage_true)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_coverage_true), 100*max(vertex_coverage_true)))

cat(sprintf("PLUGIN Coverage:\n"))
cat(sprintf("  Mean:   %.2f%% (average across all vertices)\n", 100*mean(vertex_coverage_plugin)))
cat(sprintf("  Median: %.2f%%\n", 100*median(vertex_coverage_plugin)))
cat(sprintf("  SD:     %.2f%%\n", 100*sd(vertex_coverage_plugin)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100*min(vertex_coverage_plugin), 100*max(vertex_coverage_plugin)))

# Show some example vertices
cat("--- Example Vertices ---\n")
example_vertices <- c(1, 50, 100, 250, 400, 450, 500)
cat(sprintf("%-8s %-15s %-15s %-10s\n", "Vertex", "Cov(TRUE)%", "Cov(PLUGIN)%", "NAs"))
cat(paste(rep("-", 55), collapse=""), "\n")
for (v in example_vertices) {
  cat(sprintf("%-8d %-15.1f %-15.1f %-10d\n",
              v, 100*vertex_coverage_true[v], 100*vertex_coverage_plugin[v],
              vertex_n_na_true[v]))
}
cat("\n")

# Identify vertices with extreme coverage
cat("--- Vertices with Lowest Coverage (TRUE) ---\n")
worst_vertices <- order(vertex_coverage_true)[1:10]
cat(sprintf("%-8s %-15s %-10s\n", "Vertex", "Coverage%", "NAs"))
cat(paste(rep("-", 40), collapse=""), "\n")
for (v in worst_vertices) {
  cat(sprintf("%-8d %-15.1f %-10d\n", v, 100*vertex_coverage_true[v], vertex_n_na_true[v]))
}
cat("\n")

cat("--- Vertices with Highest Coverage (TRUE) ---\n")
best_vertices <- order(vertex_coverage_true, decreasing = TRUE)[1:10]
cat(sprintf("%-8s %-15s %-10s\n", "Vertex", "Coverage%", "NAs"))
cat(paste(rep("-", 40), collapse=""), "\n")
for (v in best_vertices) {
  cat(sprintf("%-8d %-15.1f %-10d\n", v, 100*vertex_coverage_true[v], vertex_n_na_true[v]))
}
cat("\n")

cat("============================================================================\n")
cat("  PER-REPLICATION SUMMARY STATISTICS\n")
cat("============================================================================\n\n")

cat(sprintf("SSE:             Mean=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n",
            mean(rep_SSEs), sd(rep_SSEs), min(rep_SSEs), max(rep_SSEs)))
cat(sprintf("Fit Time:        Mean=%.1fs, SD=%.1fs\n",
            mean(rep_fit_times), sd(rep_fit_times)))
cat(sprintf("Coverage (TRUE):   Mean=%.1f%%, SD=%.1f%%\n",
            100*mean(rep_coverage_true), 100*sd(rep_coverage_true)))
cat(sprintf("Coverage (PLUGIN): Mean=%.1f%%, SD=%.1f%%\n\n",
            100*mean(rep_coverage_plugin), 100*sd(rep_coverage_plugin)))

# Save aggregated results
vertex_wise_results <- list(
  n = n,
  p_cov = p_cov,
  d = d,
  n_reps = n_reps,

  # Vertex-wise coverage (length n = 500)
  vertex_coverage_true = vertex_coverage_true,
  vertex_coverage_plugin = vertex_coverage_plugin,
  vertex_n_na_true = vertex_n_na_true,
  vertex_n_na_plugin = vertex_n_na_plugin,

  # Per-vertex indicators (n × n_reps matrix)
  in_ellipse_true_matrix = in_ellipse_true_matrix,
  in_ellipse_plugin_matrix = in_ellipse_plugin_matrix,

  # Per-replication statistics
  rep_ids = rep_ids,
  rep_SSEs = rep_SSEs,
  rep_fit_times = rep_fit_times,
  rep_coverage_true = rep_coverage_true,
  rep_coverage_plugin = rep_coverage_plugin
)

output_file <- file.path(result_dir, "vertex_wise_coverage_aggregated.rds")
saveRDS(vertex_wise_results, output_file)
cat(sprintf("Saved: %s\n\n", output_file))

# Also save as CSV for easy viewing
csv_data <- data.frame(
  vertex = 1:n,
  coverage_true = vertex_coverage_true,
  coverage_plugin = vertex_coverage_plugin,
  n_na_true = vertex_n_na_true,
  n_na_plugin = vertex_n_na_plugin
)
csv_file <- file.path(result_dir, "vertex_wise_coverage_results.csv")
write.csv(csv_data, csv_file, row.names = FALSE)
cat(sprintf("Saved: %s\n\n", csv_file))

cat("============================================================================\n")
cat("  AGGREGATION COMPLETE!\n")
cat("============================================================================\n")
cat("\nTo visualize vertex-wise coverage:\n")
cat("  plot(vertex_coverage_true, type='l', ylim=c(0,1),\n")
cat("       xlab='Vertex', ylab='Coverage Rate',\n")
cat("       main='Vertex-wise Coverage Across 100 Reps')\n")
cat("  abline(h=0.95, col='red', lty=2)\n")
