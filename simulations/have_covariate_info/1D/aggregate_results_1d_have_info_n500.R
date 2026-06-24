#!/usr/bin/env Rscript
# Aggregate comprehensive 1D simulation results across 100 replications
# Computes vertex-wise coverage for Fisher vs ASE vs OSE
# WITH COVARIATE SIGNAL: Z0 = rnorm(p_cov, d)

library(ggplot2)
library(dplyr)
library(tidyr)

cat("============================================================================\n")
cat("  AGGREGATING: 1D WITH COVARIATE INFO (n=500)\n")
cat("  Z0 = rnorm (covariates contain signal), maxit=30, 16 cores\n")
cat("  Methods: FISHER-TRUE/PLUGIN, ASE-TRUE/PLUGIN, OSE-TRUE/PLUGIN\n")
cat("============================================================================\n\n")

results_dir <- "results_1d_have_info_n500"

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
p_cov <- first_result$p_cov
d <- first_result$d
tau <- first_result$tau
cat(sprintf("Number of vertices: %d\n", n))
cat(sprintf("Parameters: n=%d, p_cov=%d, d=%d, tau=%.3f\n\n", n, p_cov, d, tau))

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

# Storage for convergence
converged_vec <- logical(n_reps)
iterations_vec <- integer(n_reps)

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

  # Convergence
  converged_vec[r] <- result$converged
  iterations_vec[r] <- result$iterations
}

cat("\nComputing vertex-wise coverage rates...\n")

# --- Vertex-wise coverage rates (n vertices, each with coverage rate across 100 reps) ---
vertex_rates <- list(
  fisher_true = rowMeans(fisher_true_matrix, na.rm = TRUE),
  fisher_plugin = rowMeans(fisher_plugin_matrix, na.rm = TRUE),
  ase_true = rowMeans(ase_true_matrix, na.rm = TRUE),
  ase_plugin = rowMeans(ase_plugin_matrix, na.rm = TRUE),
  ose_true = rowMeans(ose_true_matrix, na.rm = TRUE),
  ose_plugin = rowMeans(ose_plugin_matrix, na.rm = TRUE)
)

methods <- names(vertex_rates)

# --- Overall coverage statistics (statistics ACROSS vertices) ---
cat("\n============================================================================\n")
cat("OVERALL COVERAGE STATISTICS\n")
cat("============================================================================\n\n")
cat(sprintf("n=%d, p_cov=%d, d=%d, tau=%.3f, reps=%d\n", n, p_cov, d, tau, n_reps))
cat("Scenario: Z0 = rnorm (WITH covariate information)\n\n")

for (m in methods) {
  vals <- vertex_rates[[m]]  # Coverage rates across 100 reps for each of 500 vertices
  cat(sprintf("%-20s  Mean=%5.2f%%  Median=%5.2f%%  SD=%.2f%%  [%.2f%%, %.2f%%]\n",
              m,
              100 * mean(vals, na.rm = TRUE),
              100 * median(vals, na.rm = TRUE),
              100 * sd(vals, na.rm = TRUE),
              100 * min(vals, na.rm = TRUE),
              100 * max(vals, na.rm = TRUE)))
}

cat("\nSSE Statistics:\n")
cat(sprintf("  fisher    Mean=%.4f  SD=%.4f\n", mean(sse_fisher), sd(sse_fisher)))
cat(sprintf("  ase       Mean=%.4f  SD=%.4f\n", mean(sse_ase), sd(sse_ase)))
cat(sprintf("  ose       Mean=%.4f  SD=%.4f\n", mean(sse_ose), sd(sse_ose)))

cat("\nTiming (mean ± SD per replication):\n")
cat(sprintf("  fisher fit:    %.1f ± %.1f sec\n", mean(time_fisher), sd(time_fisher)))
cat(sprintf("  ASE:           %.1f ± %.1f sec\n", mean(time_ase), sd(time_ase)))
cat(sprintf("  OSE:           %.1f ± %.1f sec\n", mean(time_ose), sd(time_ose)))
cat(sprintf("  Coverage:      %.1f ± %.1f sec\n", mean(time_coverage), sd(time_coverage)))
cat(sprintf("  Total:         %.2f ± %.2f min\n",
            mean(time_fisher + time_ase + time_ose + time_coverage)/60,
            sd(time_fisher + time_ase + time_ose + time_coverage)/60))

cat(sprintf("\nfisher convergence: %.1f%%\n\n",
            100 * mean(converged_vec)))

# Create output directory for aggregated results and plots
output_dir <- "outputs_1d_have_info_n500"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save aggregated results
aggregated <- list(
  n = n,
  p_cov = p_cov,
  d = d,
  tau = tau,
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
  ),
  convergence = list(
    converged = converged_vec,
    iterations = iterations_vec,
    convergence_rate = mean(converged_vec)
  )
)

saveRDS(aggregated, file.path(output_dir, "aggregated_1d_have_info_n500.rds"))
cat(sprintf("Aggregated results saved to: %s/aggregated_1d_have_info_n500.rds\n", output_dir))

# Save summary CSV
summary_df <- data.frame(
  Method = c("fisher_true", "fisher_plugin", "ase_true", "ase_plugin", "ose_true", "ose_plugin"),
  Mean_Coverage = 100 * c(overall_fisher_true, overall_fisher_plugin,
                          overall_ase_true, overall_ase_plugin,
                          overall_ose_true, overall_ose_plugin)
)
write.csv(summary_df, file.path(output_dir, "summary_1d_have_info_n500.csv"), row.names = FALSE)
cat(sprintf("Summary saved to: %s/summary_1d_have_info_n500.csv\n\n", output_dir))

# Save vertex-wise coverage CSV
vertex_df <- data.frame(
  vertex = 1:n,
  fisher_true = vertex_rates$fisher_true,
  fisher_plugin = vertex_rates$fisher_plugin,
  ase_true = vertex_rates$ase_true,
  ase_plugin = vertex_rates$ase_plugin,
  ose_true = vertex_rates$ose_true,
  ose_plugin = vertex_rates$ose_plugin
)
write.csv(vertex_df, file.path(output_dir, "vertex_coverage_1d_have_info_n500.csv"), row.names = FALSE)
cat(sprintf("Vertex-wise CSV saved to: %s/vertex_coverage_1d_have_info_n500.csv\n", output_dir))

# Create plots
cat("\nCreating plots...\n")

# Prepare data for plotting
methods <- c("fisher_true", "fisher_plugin", "ase_true", "ase_plugin", "ose_true", "ose_plugin")
cov_matrix <- rbind(
  fisher_true_matrix,
  fisher_plugin_matrix,
  ase_true_matrix,
  ase_plugin_matrix,
  ose_true_matrix,
  ose_plugin_matrix
)
rownames(cov_matrix) <- methods

cov_long <- as.data.frame(t(cov_matrix)) |>
  mutate(rep = 1:n_reps) |>
  pivot_longer(-rep, names_to = "Method", values_to = "Coverage") |>
  mutate(
    Coverage = Coverage * 100,
    Estimator = case_when(
      grepl("fisher", Method) ~ "Fisher",
      grepl("ase", Method) ~ "ASE",
      grepl("ose", Method) ~ "OSE"
    ),
    Precision = ifelse(grepl("true", Method), "TRUE", "PLUGIN")
  )

# Plot 1: Boxplot
p1 <- ggplot(cov_long, aes(x = Method, y = Coverage, fill = Estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c(Fisher = "#E69F00", ASE = "#55A868", OSE = "#C44E52")) +
  labs(
    title = sprintf("Coverage Rates: 1D WITH Covariate Info (n=%d, %d reps)", n, n_reps),
    subtitle = "Z0 = rnorm, B contains signal | Dashed line: 95% nominal",
    x = NULL, y = "Coverage Rate (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path(output_dir, "coverage_boxplot_1d_have_info_n500.pdf"), width = 10, height = 6)
print(p1)
dev.off()

# Plot 2: Vertex-wise coverage
vertex_plot_df <- vertex_df |>
  pivot_longer(-vertex, names_to = "Method", values_to = "Coverage") |>
  mutate(Coverage = Coverage * 100)

p2 <- ggplot(vertex_plot_df, aes(x = vertex, y = Coverage, color = Method)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c(
    fisher_true = "#E69F00", fisher_plugin = "#F0E442",
    ase_true = "#2ca02c", ase_plugin = "#98df8a",
    ose_true = "#d62728", ose_plugin = "#ff9896"
  )) +
  labs(
    title = sprintf("Vertex-wise Coverage: 1D WITH Covariate Info (n=%d, %d reps)", n, n_reps),
    subtitle = "Z0 = rnorm, B contains signal | Dashed line: 95% nominal",
    x = "Vertex Index", y = "Coverage Rate (%)", color = "Method"
  ) +
  theme_minimal(base_size = 13)

pdf(file.path(output_dir, "vertex_coverage_1d_have_info_n500.pdf"), width = 12, height = 6)
print(p2)
dev.off()

cat("Plots saved:\n")
cat(sprintf("  - %s/coverage_boxplot_1d_have_info_n500.pdf\n", output_dir))
cat(sprintf("  - %s/vertex_coverage_1d_have_info_n500.pdf\n", output_dir))

cat("\n============================================================================\n")
cat("  AGGREGATION COMPLETE\n")
cat("============================================================================\n")
