#!/usr/bin/env Rscript
# Plot vertex-wise coverage rates for n=1000 simulation
# Creates a plot with vertex indices on x-axis and coverage rates on y-axis

library(ggplot2)

cat("============================================================================\n")
cat("  PLOTTING VERTEX-WISE COVERAGE (n=1000)\n")
cat("============================================================================\n\n")

# Option 1: Read from RDS file (more complete data)
# Check multiple possible locations
rds_file <- if (file.exists("vertex_wise_results_n1000/vertex_wise_coverage_aggregated_n1000.rds")) {
  "vertex_wise_results_n1000/vertex_wise_coverage_aggregated_n1000.rds"
} else {
  "vertex_wise_coverage_aggregated_n1000.rds"
}

csv_file <- if (file.exists("vertex_wise_results_n1000/vertex_wise_coverage_results_n1000.csv")) {
  "vertex_wise_results_n1000/vertex_wise_coverage_results_n1000.csv"
} else {
  "vertex_wise_coverage_results_n1000.csv"
}

# Check which file exists
if (file.exists(rds_file)) {
  cat("Reading from RDS file...\n")
  results <- readRDS(rds_file)

  # Extract vertex-wise coverage
  plot_data <- data.frame(
    vertex = 1:results$n,
    coverage_true = results$vertex_coverage_true,
    coverage_plugin = results$vertex_coverage_plugin
  )
} else if (file.exists(csv_file)) {
  cat("Reading from CSV file...\n")
  plot_data <- read.csv(csv_file)
  names(plot_data) <- c("vertex", "coverage_true", "coverage_plugin", "n_na_true", "n_na_plugin")
} else {
  stop("Neither RDS nor CSV file found!")
}

n <- nrow(plot_data)
cat(sprintf("Loaded data for %d vertices\n\n", n))

# Print summary statistics
cat("Coverage Statistics:\n")
cat(sprintf("  TRUE:   Mean=%.2f%%, Range=[%.2f%%, %.2f%%]\n",
            100*mean(plot_data$coverage_true),
            100*min(plot_data$coverage_true),
            100*max(plot_data$coverage_true)))
cat(sprintf("  PLUGIN: Mean=%.2f%%, Range=[%.2f%%, %.2f%%]\n\n",
            100*mean(plot_data$coverage_plugin),
            100*min(plot_data$coverage_plugin),
            100*max(plot_data$coverage_plugin)))

# Create output directory for plots
plot_dir <- "plots_vertex_wise_n1000"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ============================================================================
# Plot 1: Both TRUE and PLUGIN coverage on same plot
# ============================================================================
cat("Creating combined plot...\n")

pdf(file.path(plot_dir, "vertex_wise_coverage_n1000_combined.pdf"), width = 10, height = 6)
par(mar = c(4, 4, 3, 1))

plot(plot_data$vertex, plot_data$coverage_true,
     type = "l", col = "blue", lwd = 2,
     xlab = "Vertex Index", ylab = "Coverage Rate",
     main = sprintf("Vertex-wise Coverage Rates (n=%d, 100 replications)", n),
     ylim = c(0, 1), las = 1)

lines(plot_data$vertex, plot_data$coverage_plugin,
      col = "red", lwd = 2, lty = 2)

abline(h = 0.95, col = "darkgreen", lty = 2, lwd = 2)
abline(h = mean(plot_data$coverage_true), col = "blue", lty = 3)
abline(h = mean(plot_data$coverage_plugin), col = "red", lty = 3)

legend("bottomright",
       legend = c(sprintf("TRUE G_in (mean=%.2f%%)", 100*mean(plot_data$coverage_true)),
                  sprintf("PLUGIN G_in (mean=%.2f%%)", 100*mean(plot_data$coverage_plugin)),
                  "Nominal 95%"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 2),
       lwd = 2,
       bg = "white")

dev.off()
cat(sprintf("  Saved: %s\n", file.path(plot_dir, "vertex_wise_coverage_n1000_combined.pdf")))

# ============================================================================
# Plot 2: Separate plots for TRUE and PLUGIN
# ============================================================================
cat("Creating separate plots...\n")

pdf(file.path(plot_dir, "vertex_wise_coverage_n1000_separate.pdf"), width = 10, height = 8)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# TRUE coverage
plot(plot_data$vertex, plot_data$coverage_true,
     type = "l", col = "blue", lwd = 2,
     xlab = "Vertex Index", ylab = "Coverage Rate",
     main = sprintf("TRUE G_in Coverage (n=%d)", n),
     ylim = c(0, 1), las = 1)
abline(h = 0.95, col = "darkgreen", lty = 2, lwd = 2)
abline(h = mean(plot_data$coverage_true), col = "gray40", lty = 3)
text(n*0.8, 0.05,
     sprintf("Mean = %.2f%%", 100*mean(plot_data$coverage_true)),
     cex = 1.2)

# PLUGIN coverage
plot(plot_data$vertex, plot_data$coverage_plugin,
     type = "l", col = "red", lwd = 2,
     xlab = "Vertex Index", ylab = "Coverage Rate",
     main = sprintf("PLUGIN G_in Coverage (n=%d)", n),
     ylim = c(0, 1), las = 1)
abline(h = 0.95, col = "darkgreen", lty = 2, lwd = 2)
abline(h = mean(plot_data$coverage_plugin), col = "gray40", lty = 3)
text(n*0.8, 0.05,
     sprintf("Mean = %.2f%%", 100*mean(plot_data$coverage_plugin)),
     cex = 1.2)

dev.off()
cat(sprintf("  Saved: %s\n", file.path(plot_dir, "vertex_wise_coverage_n1000_separate.pdf")))

# ============================================================================
# Plot 3: ggplot2 version (if ggplot2 is available)
# ============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("Creating ggplot2 version...\n")

  library(ggplot2)

  # Reshape data for ggplot
  plot_data_long <- data.frame(
    vertex = rep(plot_data$vertex, 2),
    coverage = c(plot_data$coverage_true, plot_data$coverage_plugin),
    method = rep(c("TRUE G_in", "PLUGIN G_in"), each = n)
  )

  p <- ggplot(plot_data_long, aes(x = vertex, y = coverage, color = method, linetype = method)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkgreen", linewidth = 1) +
    scale_color_manual(values = c("TRUE G_in" = "blue", "PLUGIN G_in" = "red")) +
    scale_linetype_manual(values = c("TRUE G_in" = "solid", "PLUGIN G_in" = "dashed")) +
    labs(
      title = sprintf("Vertex-wise Coverage Rates (n=%d, 100 replications)", n),
      x = "Vertex Index",
      y = "Coverage Rate",
      color = "Method",
      linetype = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    ylim(0, 1)

  ggsave(file.path(plot_dir, "vertex_wise_coverage_n1000_ggplot.pdf"),
         plot = p, width = 10, height = 6)
  cat(sprintf("  Saved: %s\n", file.path(plot_dir, "vertex_wise_coverage_n1000_ggplot.pdf")))
}

# ============================================================================
# Save summary statistics
# ============================================================================
cat("\nSaving summary statistics...\n")

summary_stats <- data.frame(
  statistic = c("n_vertices", "n_replications",
                "coverage_true_mean", "coverage_true_median", "coverage_true_sd",
                "coverage_true_min", "coverage_true_max",
                "coverage_plugin_mean", "coverage_plugin_median", "coverage_plugin_sd",
                "coverage_plugin_min", "coverage_plugin_max"),
  value = c(
    n, 100,
    mean(plot_data$coverage_true), median(plot_data$coverage_true), sd(plot_data$coverage_true),
    min(plot_data$coverage_true), max(plot_data$coverage_true),
    mean(plot_data$coverage_plugin), median(plot_data$coverage_plugin), sd(plot_data$coverage_plugin),
    min(plot_data$coverage_plugin), max(plot_data$coverage_plugin)
  )
)

write.csv(summary_stats, file.path(plot_dir, "coverage_summary_statistics.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(plot_dir, "coverage_summary_statistics.csv")))

cat("\n============================================================================\n")
cat("  PLOTTING COMPLETE!\n")
cat(sprintf("  All plots saved to: %s/\n", plot_dir))
cat("============================================================================\n")
