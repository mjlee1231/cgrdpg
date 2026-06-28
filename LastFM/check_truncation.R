#!/usr/bin/env Rscript
# Check how many edge probabilities are truncated outside (tau, 1-tau)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript check_truncation.R <dimension>")
}

d <- as.integer(args[1])
tau <- 0.001

cat(sprintf("Analyzing truncation for d=%d, tau=%.6f\n", d, tau))
cat(sprintf("Valid range: (%.6f, %.6f)\n\n", tau, 1-tau))

# Load results
res <- readRDS(sprintf("d%d/results/lastfm_d%d_complete_results.rds", d, d))
n <- res$data$n

cat(sprintf("Network: n=%d nodes, %d unique pairs\n\n", n, n*(n-1)/2))

# Function to compute edge probabilities and count truncations
check_method <- function(X, method_name) {
  cat(sprintf("=== %s ===\n", method_name))

  # Compute all edge probabilities: s_ij = x_i^T x_j
  P <- X %*% t(X)

  # Extract upper triangle (unique pairs, no self-loops)
  probs <- P[upper.tri(P)]

  # Count truncations
  below_tau <- sum(probs < tau)
  above_1mtau <- sum(probs > (1 - tau))
  within_range <- sum(probs >= tau & probs <= (1 - tau))

  pct_below <- 100 * below_tau / length(probs)
  pct_above <- 100 * above_1mtau / length(probs)
  pct_within <- 100 * within_range / length(probs)

  cat(sprintf("  Below %.6f:     %7d pairs (%.2f%%)\n", tau, below_tau, pct_below))
  cat(sprintf("  Above %.6f:     %7d pairs (%.2f%%)\n", 1-tau, above_1mtau, pct_above))
  cat(sprintf("  Within range:      %7d pairs (%.2f%%)\n", within_range, pct_within))
  cat(sprintf("  Total truncated:   %7d pairs (%.2f%%)\n",
              below_tau + above_1mtau, pct_below + pct_above))

  # Additional statistics
  cat(sprintf("\n  Edge probability range: [%.6f, %.6f]\n", min(probs), max(probs)))
  cat(sprintf("  Mean: %.6f, Median: %.6f\n", mean(probs), median(probs)))
  cat(sprintf("  Q1: %.6f, Q3: %.6f\n\n", quantile(probs, 0.25), quantile(probs, 0.75)))

  return(list(
    below = below_tau,
    above = above_1mtau,
    within = within_range,
    pct_truncated = pct_below + pct_above,
    min = min(probs),
    max = max(probs),
    mean = mean(probs),
    median = median(probs)
  ))
}

# Check all three methods
fisher_stats <- check_method(res$latent_positions$fisher, "MODIFIED FISHER")
ase_stats <- check_method(res$latent_positions$ase, "ASE")
ose_stats <- check_method(res$latent_positions$ose, "OSE")

# Summary comparison
cat("============================================================================\n")
cat("SUMMARY: Truncation Comparison\n")
cat("============================================================================\n\n")
cat(sprintf("%-15s %12s %12s %12s\n",
            "Method", "Below tau", "Above 1-tau", "% Truncated"))
cat(sprintf("%-15s %12d %12d %12.2f%%\n",
            "Fisher", fisher_stats$below, fisher_stats$above, fisher_stats$pct_truncated))
cat(sprintf("%-15s %12d %12d %12.2f%%\n",
            "ASE", ase_stats$below, ase_stats$above, ase_stats$pct_truncated))
cat(sprintf("%-15s %12d %12d %12.2f%%\n",
            "OSE", ose_stats$below, ose_stats$above, ose_stats$pct_truncated))
cat("\n")
