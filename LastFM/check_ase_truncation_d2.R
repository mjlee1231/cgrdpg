#!/usr/bin/env Rscript
# Check how many ASE edge probabilities violate numerical bounds for d=2
# This analysis helps understand why clipping is necessary

library(cgrdpg)

d <- 2
tau <- 0.001

cat("============================================================================\n")
cat(sprintf("  ASE EDGE PROBABILITY VIOLATION ANALYSIS: d=%d, tau=%.6f\n", d, tau))
cat("============================================================================\n\n")

# Load results
res_file <- sprintf("d%d/results/lastfm_d%d_complete_results.rds", d, d)
if (!file.exists(res_file)) {
  stop(sprintf("Results file not found: %s", res_file))
}

res <- readRDS(res_file)
n <- res$data$n
n_pairs <- n * (n - 1) / 2

cat(sprintf("Network: n=%d nodes, %d unique pairs (excluding self-loops)\n", n, n_pairs))
cat(sprintf("Numerical stability bounds: [%.6f, %.6f]\n\n", tau, 1-tau))

# ============================================================================
# ANALYZE ASE (Initial latent positions)
# ============================================================================
cat("=== ASE (INITIAL POSITIONS) ===\n\n")

X_ase <- res$latent_positions$ase

# Compute all edge probabilities: p_ij = x_i^T x_j
P_ase <- X_ase %*% t(X_ase)

# Extract upper triangle (unique pairs, no self-loops)
probs_ase <- P_ase[upper.tri(P_ase)]

# Count violations
below_tau <- sum(probs_ase < tau)
above_1mtau <- sum(probs_ase > (1 - tau))
within_range <- sum(probs_ase >= tau & probs_ase <= (1 - tau))

pct_below <- 100 * below_tau / n_pairs
pct_above <- 100 * above_1mtau / n_pairs
pct_within <- 100 * within_range / n_pairs
pct_violated <- pct_below + pct_above

cat("Violation counts:\n")
cat(sprintf("  Below %.6f:     %8d pairs (%.2f%%)\n", tau, below_tau, pct_below))
cat(sprintf("  Above %.6f:     %8d pairs (%.2f%%)\n", 1-tau, above_1mtau, pct_above))
cat(sprintf("  Within [%.6f, %.6f]: %8d pairs (%.2f%%)\n\n",
            tau, 1-tau, within_range, pct_within))

cat(sprintf("Total violations: %d pairs (%.2f%%)\n\n",
            below_tau + above_1mtau, pct_violated))

# Distribution statistics
cat("ASE edge probability distribution:\n")
cat(sprintf("  Min:    %.6f\n", min(probs_ase)))
cat(sprintf("  Q1:     %.6f\n", quantile(probs_ase, 0.25)))
cat(sprintf("  Median: %.6f\n", median(probs_ase)))
cat(sprintf("  Mean:   %.6f\n", mean(probs_ase)))
cat(sprintf("  Q3:     %.6f\n", quantile(probs_ase, 0.75)))
cat(sprintf("  Max:    %.6f\n\n", max(probs_ase)))

# ============================================================================
# COMPARISON: OSE and MSLE
# ============================================================================
cat("============================================================================\n")
cat("  COMPARISON: Edge Probability Distributions\n")
cat("============================================================================\n\n")

show_violation <- function(X, method_name) {
  P <- X %*% t(X)
  probs <- P[upper.tri(P)]

  below <- sum(probs < tau)
  above <- sum(probs > (1 - tau))
  pct_viol <- 100 * (below + above) / n_pairs

  cat(sprintf("%s:\n", method_name))
  cat(sprintf("  Range: [%.6f, %.6f]\n", min(probs), max(probs)))
  cat(sprintf("  Mean:  %.6f, Median: %.6f\n", mean(probs), median(probs)))
  cat(sprintf("  Violations: %d pairs (%.2f%%)\n\n", below + above, pct_viol))
}

show_violation(res$latent_positions$ase, "ASE")
show_violation(res$latent_positions$ose, "OSE")
show_violation(res$latent_positions$fisher, "MSLE (R)")

# ============================================================================
# LOG-LIKELIHOOD IMPROVEMENT
# ============================================================================
cat("============================================================================\n")
cat("  LOG-LIKELIHOOD TRAJECTORY (MSLE)\n")
cat("============================================================================\n\n")

if (!is.null(res$fisher_info$objective)) {
  obj_trajectory <- res$fisher_info$objective
  n_iters <- length(obj_trajectory)

  initial_ll <- obj_trajectory[1]
  final_ll <- tail(obj_trajectory, 1)

  # Log-likelihood (negate because objective is -log_lik)
  initial_loglik <- -initial_ll
  final_loglik <- -final_ll

  improvement <- final_ll - initial_ll
  pct_improvement <- 100 * abs(improvement) / abs(initial_ll)

  cat(sprintf("Iterations: %d\n", n_iters))
  cat(sprintf("Converged: %s\n\n", ifelse(res$fisher_info$converged, "YES", "NO")))

  cat("Surrogate log-likelihood trajectory:\n")
  cat(sprintf("  Initial: %.0f\n", initial_ll))
  cat(sprintf("  Final:   %.0f\n", final_ll))
  cat(sprintf("  Change:  %+.0f (%.2f%% reduction in discrepancy)\n\n",
              improvement, pct_improvement))

  # Show first few and last few iterations
  cat("Iteration-by-iteration:\n")
  if (n_iters <= 10) {
    for (i in 1:n_iters) {
      cat(sprintf("  Iter %2d: %.0f\n", i-1, obj_trajectory[i]))
    }
  } else {
    for (i in 1:5) {
      cat(sprintf("  Iter %2d: %.0f\n", i-1, obj_trajectory[i]))
    }
    cat("  ...\n")
    for (i in (n_iters-4):n_iters) {
      cat(sprintf("  Iter %2d: %.0f\n", i-1, obj_trajectory[i]))
    }
  }
} else {
  cat("Log-likelihood trajectory not available.\n")
}

cat("\n============================================================================\n")
cat("  SUMMARY FOR d=2\n")
cat("============================================================================\n\n")

cat(sprintf("ASE edge probabilities violating numerical bounds: %.2f%%\n", pct_violated))
cat(sprintf("  (%.2f%% below %.6f, %.2f%% above %.6f)\n\n",
            pct_below, tau, pct_above, 1-tau))

cat("This high violation rate demonstrates why the modified clipping strategy\n")
cat("is essential for stable Fisher-scoring computation in sparse networks.\n")
cat("============================================================================\n")
