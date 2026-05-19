#!/usr/bin/env Rscript
# Aggregate vertex-wise coverage results for 3D GRDPG: n=500

library(ggplot2)
library(dplyr)

cat("============================================================================\n")
cat("  AGGREGATING VERTEX-WISE COVERAGE RESULTS: 3D GRDPG, n=500\n")
cat("============================================================================\n\n")

# Load all results
results_dir <- "results_3d_n500"
if (!dir.exists(results_dir)) {
  stop(sprintf("Results directory '%s' not found!", results_dir))
}

files <- list.files(results_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
n_reps <- length(files)

cat(sprintf("Found %d replication files\n\n", n_reps))

if (n_reps == 0) {
  stop("No result files found!")
}

# Load all results
all_results <- lapply(files, readRDS)

# Extract summary statistics
n <- all_results[[1]]$n
d <- all_results[[1]]$d
p_cov <- all_results[[1]]$p_cov
tau <- all_results[[1]]$tau

# Aggregate coverage rates
coverage_true_all <- sapply(all_results, function(x) x$coverage_true)
coverage_plugin_all <- sapply(all_results, function(x) x$coverage_plugin)

# Aggregate vertex-wise coverage
vertex_coverage_true <- matrix(0, n, n_reps)
vertex_coverage_plugin <- matrix(0, n, n_reps)

for (i in 1:n_reps) {
  vertex_coverage_true[, i] <- all_results[[i]]$in_ellipse_true
  vertex_coverage_plugin[, i] <- all_results[[i]]$in_ellipse_plugin
}

# Compute vertex-wise coverage rates (across replications)
vertex_rate_true <- rowMeans(vertex_coverage_true, na.rm = TRUE)
vertex_rate_plugin <- rowMeans(vertex_coverage_plugin, na.rm = TRUE)

# Overall statistics
cat("============================================================================\n")
cat("OVERALL COVERAGE STATISTICS (across replications)\n")
cat("============================================================================\n\n")

cat(sprintf("Number of replications: %d\n", n_reps))
cat(sprintf("n=%d, p_cov=%d, d=%d, tau=%.3f\n\n", n, p_cov, d, tau))

cat("Coverage rates (TRUE G_in):\n")
cat(sprintf("  Mean:   %.2f%%\n", 100 * mean(coverage_true_all)))
cat(sprintf("  Median: %.2f%%\n", 100 * median(coverage_true_all)))
cat(sprintf("  SD:     %.2f%%\n", 100 * sd(coverage_true_all)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100 * min(coverage_true_all), 100 * max(coverage_true_all)))

cat("Coverage rates (PLUGIN G_in):\n")
cat(sprintf("  Mean:   %.2f%%\n", 100 * mean(coverage_plugin_all)))
cat(sprintf("  Median: %.2f%%\n", 100 * median(coverage_plugin_all)))
cat(sprintf("  SD:     %.2f%%\n", 100 * sd(coverage_plugin_all)))
cat(sprintf("  Range:  [%.2f%%, %.2f%%]\n\n",
            100 * min(coverage_plugin_all), 100 * max(coverage_plugin_all)))

# Timing statistics
fit_times <- sapply(all_results, function(x) x$fit_time)
coverage_times <- sapply(all_results, function(x) x$coverage_time)

cat("Timing statistics:\n")
cat(sprintf("  Mean fit time:      %.1f sec (%.2f min)\n", mean(fit_times), mean(fit_times)/60))
cat(sprintf("  Mean coverage time: %.1f sec (%.2f min)\n", mean(coverage_times), mean(coverage_times)/60))
cat(sprintf("  Total mean time:    %.2f min per replication\n\n", mean(fit_times + coverage_times)/60))

# SSE statistics
sse_all <- sapply(all_results, function(x) x$SSE)
cat("SSE statistics:\n")
cat(sprintf("  Mean:   %.4f\n", mean(sse_all)))
cat(sprintf("  Median: %.4f\n", median(sse_all)))
cat(sprintf("  SD:     %.4f\n\n", sd(sse_all)))

# Convergence statistics
converged_all <- sapply(all_results, function(x) x$converged)
cat(sprintf("Convergence: %d/%d (%.1f%%)\n\n", sum(converged_all), n_reps, 100*mean(converged_all)))

# Save aggregated results
aggregated <- list(
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
  n_reps = n_reps,
  coverage_true_all = coverage_true_all,
  coverage_plugin_all = coverage_plugin_all,
  vertex_rate_true = vertex_rate_true,
  vertex_rate_plugin = vertex_rate_plugin,
  sse_all = sse_all,
  fit_times = fit_times,
  coverage_times = coverage_times,
  converged_all = converged_all
)

saveRDS(aggregated, "aggregated_results_3d_n500.rds")
cat("Aggregated results saved to: aggregated_results_3d_n500.rds\n\n")

# Create summary CSV
summary_df <- data.frame(
  Method = c("TRUE", "PLUGIN"),
  Mean_Coverage = c(mean(coverage_true_all), mean(coverage_plugin_all)),
  SD_Coverage = c(sd(coverage_true_all), sd(coverage_plugin_all)),
  Min_Coverage = c(min(coverage_true_all), min(coverage_plugin_all)),
  Max_Coverage = c(max(coverage_true_all), max(coverage_plugin_all))
)

write.csv(summary_df, "summary_3d_n500.csv", row.names = FALSE)
cat("Summary saved to: summary_3d_n500.csv\n\n")

# Plot coverage distributions
cat("Creating plots...\n")

# 1. Boxplot of coverage rates
pdf("coverage_boxplot_3d_n500.pdf", width = 8, height = 6)
boxplot_data <- data.frame(
  Coverage = c(coverage_true_all * 100, coverage_plugin_all * 100),
  Method = rep(c("TRUE", "PLUGIN"), each = n_reps)
)
p <- ggplot(boxplot_data, aes(x = Method, y = Coverage, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  labs(title = sprintf("Coverage Rates: 3D GRDPG (n=%d, %d reps)", n, n_reps),
       y = "Coverage Rate (%)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
dev.off()

# 2. Vertex-wise coverage plot
pdf("vertex_coverage_3d_n500.pdf", width = 10, height = 6)
vertex_df <- data.frame(
  Vertex = 1:n,
  TRUE_Coverage = vertex_rate_true * 100,
  PLUGIN_Coverage = vertex_rate_plugin * 100
)
p <- ggplot(vertex_df) +
  geom_line(aes(x = Vertex, y = TRUE_Coverage, color = "TRUE"), linewidth = 0.8) +
  geom_line(aes(x = Vertex, y = PLUGIN_Coverage, color = "PLUGIN"), linewidth = 0.8) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  labs(title = sprintf("Vertex-wise Coverage Rates: 3D GRDPG (n=%d, %d reps)", n, n_reps),
       x = "Vertex Index", y = "Coverage Rate (%)", color = "Method") +
  theme_minimal()
print(p)
dev.off()

cat("Plots saved:\n")
cat("  - coverage_boxplot_3d_n500.pdf\n")
cat("  - vertex_coverage_3d_n500.pdf\n")

cat("\n============================================================================\n")
cat("AGGREGATION COMPLETE\n")
cat("============================================================================\n")
