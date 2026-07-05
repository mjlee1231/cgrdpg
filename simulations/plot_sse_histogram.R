#!/usr/bin/env Rscript
# Plot histogram of SSEs from simulation results

library(ggplot2)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_sse_histogram.R <results_file.rds> [output_prefix]\n")
  cat("Example: Rscript plot_sse_histogram.R results_3d_ase_ose_cgrdpg_n1000.rds sse_n1000\n")
  quit(status = 1)
}

results_file <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "sse_histogram"

# Check file exists
if (!file.exists(results_file)) {
  cat(sprintf("Error: File '%s' not found\n", results_file))
  quit(status = 1)
}

cat(sprintf("Loading results from: %s\n", results_file))
results <- readRDS(results_file)

cat(sprintf("Results structure: %s\n", paste(names(results), collapse=", ")))

# Extract SSE data
# Assuming structure: results$sse_ase, results$sse_ose, results$sse_cgrdpg
# Or results might be a list of replications

# Try to infer structure
if ("sse_ase" %in% names(results)) {
  # Structure: direct SSE fields
  sse_data <- data.frame(
    SSE = c(results$sse_ase, results$sse_ose, results$sse_cgrdpg),
    Method = rep(c("ASE", "OSE", "CGRDPG"),
                 each = length(results$sse_ase))
  )
} else if (is.list(results) && "sse_ase" %in% names(results[[1]])) {
  # Structure: list of replications
  n_reps <- length(results)
  sse_ase <- sapply(results, function(x) x$sse_ase)
  sse_ose <- sapply(results, function(x) x$sse_ose)
  sse_cgrdpg <- sapply(results, function(x) x$sse_cgrdpg)

  sse_data <- data.frame(
    SSE = c(sse_ase, sse_ose, sse_cgrdpg),
    Method = rep(c("ASE", "OSE", "CGRDPG"), each = n_reps)
  )
} else {
  cat("Error: Cannot identify SSE data structure\n")
  cat("Available fields:\n")
  print(str(results, max.level = 1))
  quit(status = 1)
}

cat(sprintf("\nLoaded %d SSE values across %d methods\n",
            nrow(sse_data), length(unique(sse_data$Method))))

# Summary statistics
cat("\nSummary statistics by method:\n")
print(aggregate(SSE ~ Method, data = sse_data,
                FUN = function(x) c(
                  n = length(x),
                  mean = mean(x),
                  sd = sd(x),
                  median = median(x),
                  min = min(x),
                  max = max(x)
                )))

# Create histogram
cat("\nCreating histogram...\n")

# Use Cairo for HPC compatibility
options(bitmapType = 'cairo')

p <- ggplot(sse_data, aes(x = SSE, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of Sum of Squared Errors (SSE)",
    subtitle = sprintf("Based on %s", basename(results_file)),
    x = "SSE",
    y = "Frequency"
  ) +
  theme(legend.position = "none")

# Save plot
output_file <- sprintf("%s.png", output_prefix)
ggsave(output_file, p, width = 8, height = 10, dpi = 300)
cat(sprintf("✓ Histogram saved to: %s\n", output_file))

# Also create overlaid histogram
p2 <- ggplot(sse_data, aes(x = SSE, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() +
  labs(
    title = "Distribution of Sum of Squared Errors (SSE) - Overlaid",
    subtitle = sprintf("Based on %s", basename(results_file)),
    x = "SSE",
    y = "Frequency"
  ) +
  scale_fill_brewer(palette = "Set2")

output_file2 <- sprintf("%s_overlaid.png", output_prefix)
ggsave(output_file2, p2, width = 10, height = 6, dpi = 300)
cat(sprintf("✓ Overlaid histogram saved to: %s\n", output_file2))

cat("\nDone!\n")
