#!/usr/bin/env Rscript
# Aggregate vertex-wise coverage results: 3D GRDPG, ASE vs OSE vs cgrdpg, n=1000

library(ggplot2)
library(dplyr)
library(tidyr)

cat("============================================================================\n")
cat("  AGGREGATING: 3D GRDPG ASE vs OSE vs cgrdpg, n=1000\n")
cat("  Methods: cgrdpg-TRUE, cgrdpg-PLUGIN, ASE-TRUE, ASE-PLUGIN, OSE-TRUE, OSE-PLUGIN\n")
cat("============================================================================\n\n")

results_dir <- "results_3d_ase_ose_n1000"
if (!dir.exists(results_dir)) stop(sprintf("Directory '%s' not found!", results_dir))

files  <- list.files(results_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
n_reps <- length(files)
cat(sprintf("Found %d replication files\n\n", n_reps))
if (n_reps == 0) stop("No result files found!")

all_results <- lapply(files, readRDS)

# --- Extract meta ---
n     <- all_results[[1]]$n
p_cov <- all_results[[1]]$p_cov
d     <- all_results[[1]]$d
tau   <- all_results[[1]]$tau
methods <- colnames(all_results[[1]]$results_mat)

# --- Overall coverage per rep (6 x n_reps matrix) ---
cov_matrix <- sapply(all_results, function(x) x$overall_cov)  # 6 x n_reps
rownames(cov_matrix) <- methods

cat("============================================================================\n")
cat("OVERALL COVERAGE STATISTICS\n")
cat("============================================================================\n\n")
cat(sprintf("n=%d, p_cov=%d, d=%d, tau=%.3f, reps=%d\n\n", n, p_cov, d, tau, n_reps))

for (m in methods) {
  vals <- cov_matrix[m, ]
  cat(sprintf("%-20s  Mean=%5.2f%%  Median=%5.2f%%  SD=%.2f%%  [%.2f%%, %.2f%%]\n",
              m,
              100 * mean(vals),
              100 * median(vals),
              100 * sd(vals),
              100 * min(vals),
              100 * max(vals)))
}

# --- Vertex-wise coverage rates (n x n_reps per method) ---
cat("\nComputing vertex-wise coverage rates...\n")
vertex_rates <- list()
for (m in methods) {
  mat <- sapply(all_results, function(x) x$results_mat[, m])  # n x n_reps
  vertex_rates[[m]] <- rowMeans(mat, na.rm = TRUE)
}

# --- SSE ---
cat("\nSSE Statistics:\n")
for (est in c("cgrdpg", "ase", "ose")) {
  sse_vals <- sapply(all_results, function(x) x$sse[est])
  cat(sprintf("  %-8s  Mean=%.4f  SD=%.4f\n", est, mean(sse_vals), sd(sse_vals)))
}

# --- Timing ---
cat("\nTiming (mean per replication):\n")
cgrdpg_t  <- mean(sapply(all_results, function(x) x$timing$cgrdpg_time))
ase_t     <- mean(sapply(all_results, function(x) x$timing$ase_time))
ose_t     <- mean(sapply(all_results, function(x) x$timing$ose_time))
cov_t     <- mean(sapply(all_results, function(x) x$timing$coverage_time))
total_t   <- mean(sapply(all_results, function(x) x$timing$rep_time_min))
cat(sprintf("  cgrdpg fit:    %.1f sec\n", cgrdpg_t))
cat(sprintf("  ASE:           %.1f sec\n", ase_t))
cat(sprintf("  OSE:           %.1f sec\n", ose_t))
cat(sprintf("  Coverage:      %.1f sec\n", cov_t))
cat(sprintf("  Total:         %.2f min\n", total_t))

# --- Convergence ---
conv <- mean(sapply(all_results, function(x) x$converged))
cat(sprintf("\ncgrdpg convergence: %.1f%%\n\n", 100 * conv))

# --- Save aggregated results ---
aggregated <- list(
  n            = n,
  p_cov        = p_cov,
  d            = d,
  tau          = tau,
  n_reps       = n_reps,
  methods      = methods,
  cov_matrix   = cov_matrix,    # 6 x n_reps
  vertex_rates = vertex_rates,  # list of 6 vectors length n
  sse_all      = sapply(all_results, function(x) x$sse),
  convergence  = conv
)
saveRDS(aggregated, "aggregated_3d_ase_ose_n1000.rds")
cat("Aggregated results saved to: aggregated_3d_ase_ose_n1000.rds\n")

# --- Summary CSV ---
summary_df <- data.frame(
  Method         = methods,
  Mean_Coverage  = 100 * rowMeans(cov_matrix),
  SD_Coverage    = 100 * apply(cov_matrix, 1, sd),
  Min_Coverage   = 100 * apply(cov_matrix, 1, min),
  Max_Coverage   = 100 * apply(cov_matrix, 1, max)
)
write.csv(summary_df, "summary_3d_ase_ose_n1000.csv", row.names = FALSE)
cat("Summary saved to: summary_3d_ase_ose_n1000.csv\n\n")

# --- Plot 1: Boxplot of overall coverage per method ---
cat("Creating plots...\n")

cov_long <- as.data.frame(t(cov_matrix)) |>
  mutate(rep = 1:n_reps) |>
  pivot_longer(-rep, names_to = "Method", values_to = "Coverage") |>
  mutate(
    Coverage = Coverage * 100,
    Estimator = case_when(
      grepl("cgrdpg", Method) ~ "cgrdpg",
      grepl("ase",    Method) ~ "ASE",
      grepl("ose",    Method) ~ "OSE"
    ),
    Precision = ifelse(grepl("true", Method), "TRUE", "PLUGIN")
  )

p1 <- ggplot(cov_long, aes(x = Method, y = Coverage, fill = Estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c(cgrdpg = "#4C72B0", ASE = "#55A868", OSE = "#C44E52")) +
  labs(
    title = sprintf("Coverage Rates: 3D GRDPG (n=%d, %d reps)", n, n_reps),
    subtitle = "S = diag(1,1,-1), p_cov=500 | Dashed line: 95% nominal",
    x = NULL, y = "Coverage Rate (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf("coverage_boxplot_3d_ase_ose_n1000.pdf", width = 10, height = 6)
print(p1)
dev.off()

# --- Plot 2: Vertex-wise coverage ---
vertex_df <- data.frame(
  Vertex       = 1:n,
  cgrdpg_true   = vertex_rates[["cgrdpg_true"]]   * 100,
  cgrdpg_plugin = vertex_rates[["cgrdpg_plugin"]] * 100,
  ase_true      = vertex_rates[["ase_true"]]       * 100,
  ase_plugin    = vertex_rates[["ase_plugin"]]     * 100,
  ose_true      = vertex_rates[["ose_true"]]       * 100,
  ose_plugin    = vertex_rates[["ose_plugin"]]     * 100
) |>
  pivot_longer(-Vertex, names_to = "Method", values_to = "Coverage")

p2 <- ggplot(vertex_df, aes(x = Vertex, y = Coverage, color = Method)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c(
    cgrdpg_true   = "#1f77b4", cgrdpg_plugin = "#aec7e8",
    ase_true      = "#2ca02c", ase_plugin    = "#98df8a",
    ose_true      = "#d62728", ose_plugin    = "#ff9896"
  )) +
  labs(
    title   = sprintf("Vertex-wise Coverage: 3D GRDPG (n=%d, %d reps)", n, n_reps),
    subtitle = "S = diag(1,1,-1), p_cov=500 | Dashed line: 95% nominal",
    x = "Vertex Index", y = "Coverage Rate (%)", color = "Method"
  ) +
  theme_minimal(base_size = 13)

pdf("vertex_coverage_3d_ase_ose_n1000.pdf", width = 12, height = 6)
print(p2)
dev.off()

# --- Plot 3: TRUE vs PLUGIN side-by-side per estimator ---
p3 <- ggplot(cov_long, aes(x = Precision, y = Coverage, fill = Precision)) +
  geom_boxplot() +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  facet_wrap(~ Estimator, nrow = 1) +
  scale_fill_manual(values = c(TRUE = "#4C72B0", PLUGIN = "#DD8452")) +
  labs(
    title   = sprintf("TRUE vs PLUGIN Coverage by Estimator: 3D GRDPG (n=%d, %d reps)", n, n_reps),
    x = NULL, y = "Coverage Rate (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

pdf("true_vs_plugin_3d_ase_ose_n1000.pdf", width = 10, height = 5)
print(p3)
dev.off()

cat("Plots saved:\n")
cat("  - coverage_boxplot_3d_ase_ose_n1000.pdf\n")
cat("  - vertex_coverage_3d_ase_ose_n1000.pdf\n")
cat("  - true_vs_plugin_3d_ase_ose_n1000.pdf\n")
cat("\n============================================================================\n")
cat("AGGREGATION COMPLETE\n")
cat("============================================================================\n")
