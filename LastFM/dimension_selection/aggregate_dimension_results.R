#!/usr/bin/env Rscript
# Aggregate dimension selection results from parallel jobs

# Use Cairo for graphics on HPC (no X11 display)
options(bitmapType = 'cairo')

library(ggplot2)

cat("============================================================================\n")
cat("  AGGREGATING DIMENSION SELECTION RESULTS\n")
cat("============================================================================\n\n")

# Check for result files (tau=0.001)
result_files <- paste0("results_tau0.001/results_d", 1:5, ".rds")
found_files <- file.exists(result_files)

if (!all(found_files)) {
  cat("ERROR: Not all result files found!\n")
  cat("Missing files:\n")
  for (i in which(!found_files)) {
    cat(sprintf("  - results_tau0.001/results_d%d.rds\n", i))
  }
  stop("Cannot aggregate incomplete results")
}

cat("Found all 5 result files\n\n")

# Load all results (tau=0.001)
all_results <- lapply(1:5, function(d) readRDS(sprintf("results_tau0.001/results_d%d.rds", d)))

# Extract summary table
results <- data.frame(
  d = sapply(all_results, function(x) x$d),
  neg_log_lik = sapply(all_results, function(x) x$neg_log_lik),
  n_params = sapply(all_results, function(x) x$n_params),
  BIC = sapply(all_results, function(x) x$BIC),
  AIC = sapply(all_results, function(x) x$AIC),
  converged = sapply(all_results, function(x) x$converged),
  iters = sapply(all_results, function(x) x$iters),
  time_sec = sapply(all_results, function(x) x$time_sec)
)

# ============================================================================
# DISPLAY RESULTS
# ============================================================================
cat("============================================================================\n")
cat("DIMENSION SELECTION RESULTS\n")
cat("============================================================================\n\n")

print(results, row.names = FALSE)

# Find optimal dimensions
optimal_bic <- results$d[which.min(results$BIC)]
optimal_aic <- results$d[which.min(results$AIC)]

cat(sprintf("\n\nOPTIMAL DIMENSION:\n"))
cat(sprintf("  By BIC: d = %d (BIC = %.2f)\n", optimal_bic, min(results$BIC, na.rm=TRUE)))
cat(sprintf("  By AIC: d = %d (AIC = %.2f)\n", optimal_aic, min(results$AIC, na.rm=TRUE)))

# Save aggregated results (tau=0.001)
saveRDS(list(
  results = results,
  optimal_bic = optimal_bic,
  optimal_aic = optimal_aic,
  all_fits = all_results
), "results_tau0.001/dimension_selection_aggregated.rds")

write.csv(results, "results_tau0.001/dimension_selection_summary.csv", row.names = FALSE)

cat("\n\nResults saved to:\n")
cat("  - results_tau0.001/dimension_selection_aggregated.rds (complete results)\n")
cat("  - results_tau0.001/dimension_selection_summary.csv (summary table)\n\n")

# ============================================================================
# PLOTS
# ============================================================================
cat("Creating plots...\n")

# Plot 1: BIC and AIC vs dimension
results_long <- data.frame(
  d = rep(results$d, 2),
  value = c(results$BIC, results$AIC),
  criterion = rep(c("BIC", "AIC"), each = nrow(results))
)

p1 <- ggplot(results_long, aes(x = d, y = value, color = criterion, group = criterion)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = optimal_bic, linetype = "dashed", color = "#F8766D", alpha = 0.7) +
  geom_vline(xintercept = optimal_aic, linetype = "dashed", color = "#00BFC4", alpha = 0.7) +
  scale_x_continuous(breaks = 1:5) +
  scale_color_manual(values = c("BIC" = "#F8766D", "AIC" = "#00BFC4")) +
  labs(
    title = "Model Selection: BIC and AIC vs Embedding Dimension",
    subtitle = sprintf("Optimal: BIC selects d=%d, AIC selects d=%d", optimal_bic, optimal_aic),
    x = "Embedding Dimension (d)",
    y = "Information Criterion",
    color = "Criterion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave("results_tau0.001/dimension_selection_criteria.pdf", p1, width = 8, height = 6)
ggsave("results_tau0.001/dimension_selection_criteria.png", p1, width = 8, height = 6, dpi = 300)

# Plot 2: Negative log-likelihood vs dimension
p2 <- ggplot(results, aes(x = d, y = neg_log_lik)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = 1:5) +
  labs(
    title = "Negative Log Pseudo-Likelihood vs Dimension",
    subtitle = "Lower is better (better fit to data)",
    x = "Embedding Dimension (d)",
    y = "Negative Log Pseudo-Likelihood"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("results_tau0.001/dimension_selection_loglik.pdf", p2, width = 8, height = 6)
ggsave("results_tau0.001/dimension_selection_loglik.png", p2, width = 8, height = 6, dpi = 300)

# Plot 3: Number of parameters vs dimension
p3 <- ggplot(results, aes(x = d, y = n_params)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  scale_x_continuous(breaks = 1:5) +
  labs(
    title = "Model Complexity vs Dimension",
    subtitle = sprintf("Parameters = n*d + p*d - d²"),
    x = "Embedding Dimension (d)",
    y = "Number of Parameters"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("results_tau0.001/dimension_selection_params.pdf", p3, width = 8, height = 6)
ggsave("results_tau0.001/dimension_selection_params.png", p3, width = 8, height = 6, dpi = 300)

# Plot 4: Computation time vs dimension
p4 <- ggplot(results, aes(x = d, y = time_sec / 3600)) +
  geom_line(linewidth = 1, color = "purple") +
  geom_point(size = 3, color = "purple") +
  scale_x_continuous(breaks = 1:5) +
  labs(
    title = "Computation Time vs Dimension",
    subtitle = "32 cores parallel per dimension (tau=0.001)",
    x = "Embedding Dimension (d)",
    y = "Time (hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("results_tau0.001/dimension_selection_timing.pdf", p4, width = 8, height = 6)
ggsave("results_tau0.001/dimension_selection_timing.png", p4, width = 8, height = 6, dpi = 300)

cat("\nPlots saved:\n")
cat("  - results_tau0.001/dimension_selection_criteria.pdf/png\n")
cat("  - results_tau0.001/dimension_selection_loglik.pdf/png\n")
cat("  - results_tau0.001/dimension_selection_params.pdf/png\n")
cat("  - results_tau0.001/dimension_selection_timing.pdf/png\n")

cat("\n============================================================================\n")
cat("AGGREGATION COMPLETE\n")
cat("============================================================================\n")
cat(sprintf("\nRECOMMENDATION: Use d = %d for subsequent analysis (selected by BIC)\n", optimal_bic))
cat("============================================================================\n")
