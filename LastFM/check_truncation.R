#!/usr/bin/env Rscript
# Check how many edge probabilities are truncated in Fisher method
# Note: Only Fisher uses tau for truncation. ASE/OSE are spectral methods without boundaries.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript check_truncation.R <dimension>")
}

d <- as.integer(args[1])
tau <- 0.001

cat("============================================================================\n")
cat(sprintf("  FISHER TRUNCATION ANALYSIS: d=%d, tau=%.6f\n", d, tau))
cat("============================================================================\n\n")

# Load results
res <- readRDS(sprintf("d%d/results/lastfm_d%d_complete_results.rds", d, d))
n <- res$data$n
n_pairs <- n * (n - 1) / 2

cat(sprintf("Network: n=%d nodes, %d unique pairs (excluding self-loops)\n", n, n_pairs))
cat(sprintf("Truncation bounds: (%.6f, %.6f)\n\n", tau, 1-tau))

# Analyze Fisher method
cat("=== MODIFIED FISHER METHOD ===\n\n")

X_fisher <- res$latent_positions$fisher

# Compute all edge probabilities: p_ij = x_i^T x_j
P_fisher <- X_fisher %*% t(X_fisher)

# Extract upper triangle (unique pairs, no self-loops)
probs_fisher <- P_fisher[upper.tri(P_fisher)]

# Count truncations
below_tau <- sum(probs_fisher < tau)
above_1mtau <- sum(probs_fisher > (1 - tau))
within_range <- sum(probs_fisher >= tau & probs_fisher <= (1 - tau))

pct_below <- 100 * below_tau / n_pairs
pct_above <- 100 * above_1mtau / n_pairs
pct_within <- 100 * within_range / n_pairs
pct_truncated <- pct_below + pct_above

cat("Truncation counts:\n")
cat(sprintf("  Below %.6f:     %8d pairs (%.2f%%)\n", tau, below_tau, pct_below))
cat(sprintf("  Above %.6f:     %8d pairs (%.2f%%)\n", 1-tau, above_1mtau, pct_above))
cat(sprintf("  Within [%.6f, %.6f]: %8d pairs (%.2f%%)\n\n",
            tau, 1-tau, within_range, pct_within))

cat(sprintf("Total truncated: %d pairs (%.2f%%)\n\n",
            below_tau + above_1mtau, pct_truncated))

# Distribution statistics
cat("Edge probability distribution:\n")
cat(sprintf("  Min:    %.6f\n", min(probs_fisher)))
cat(sprintf("  Q1:     %.6f\n", quantile(probs_fisher, 0.25)))
cat(sprintf("  Median: %.6f\n", median(probs_fisher)))
cat(sprintf("  Mean:   %.6f\n", mean(probs_fisher)))
cat(sprintf("  Q3:     %.6f\n", quantile(probs_fisher, 0.75)))
cat(sprintf("  Max:    %.6f\n\n", max(probs_fisher)))

# Additional context: show ASE/OSE distributions for comparison
cat("============================================================================\n")
cat("  COMPARISON: Edge Probability Distributions (for reference)\n")
cat("  Note: ASE and OSE don't use tau - these are just their natural ranges\n")
cat("============================================================================\n\n")

show_distribution <- function(X, method_name) {
  P <- X %*% t(X)
  probs <- P[upper.tri(P)]

  cat(sprintf("%s:\n", method_name))
  cat(sprintf("  Range: [%.6f, %.6f]\n", min(probs), max(probs)))
  cat(sprintf("  Mean:  %.6f, Median: %.6f\n\n", mean(probs), median(probs)))
}

show_distribution(res$latent_positions$fisher, "Modified Fisher")
show_distribution(res$latent_positions$ase, "ASE")
show_distribution(res$latent_positions$ose, "OSE")

cat("============================================================================\n")
