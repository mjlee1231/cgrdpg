#!/usr/bin/env Rscript
# Diagnostic: Check if Fisher-scoring is feasible for LastFM data
# Analyzes ASE initialization to predict Fisher-scoring performance

library(cgrdpg)
library(Matrix)

cat("============================================================================\n")
cat("  FISHER-SCORING FEASIBILITY DIAGNOSTIC\n")
cat("  Checking edge probability distribution from ASE initialization\n")
cat("============================================================================\n\n")

# ============================================================================
# PREPROCESSING (same as main analysis)
# ============================================================================
cat("Loading and preprocessing data...\n")

load("CovariatesMatrix.Rdata")
X <- as.matrix(sparseX)
edgelist <- read.csv("lastfm_asia_edges.csv") + 1
A <- matrix(0, nrow(X), nrow(X))
for (i in 1:nrow(edgelist)) A[edgelist[i,1], edgelist[i,2]] <- 1
A <- A + t(A)

label <- read.csv('lastfm_asia_target.csv')$target + 1
ind5 <- which(label != 5)
A <- A[ind5, ind5]; X <- X[ind5,]; label <- label[ind5]
X <- X[, colSums(X) > 0]

labelnew <- label; label[labelnew > 5] <- label[labelnew > 5] - 1

sizes <- summary(as.factor(label))
class_select <- names(sizes)[sizes > 300 & sizes < 1000]
ind.select <- which(label %in% class_select)
Aselect <- A[ind.select, ind.select]
Xselect <- X[ind.select,]
labelselect <- label[ind.select]

dartist_all <- colSums(X)
dartist <- colSums(Xselect)
prop <- dartist / dartist_all
prob <- 1 - round(min(nrow(Aselect)/2, 600)/ncol(Xselect), 4)
Xselect <- Xselect[, prop > quantile(prop, probs=prob, na.rm=TRUE)]

dnode <- rowSums(Xselect); like <- which(dnode > 0)
Xselect <- Xselect[like,]; Aselect <- Aselect[like,like]; labelselect <- labelselect[like]

dselect <- rowSums(Aselect)
while(sum(dselect <= 1) > 0) {
  keep <- which(dselect > 1)
  Aselect <- Aselect[keep,keep]; Xselect <- Xselect[keep,]; labelselect <- labelselect[keep]
  dselect <- rowSums(Aselect)
}

n <- nrow(Aselect)
p <- ncol(Xselect)
B <- t(Xselect)
d <- 2

cat(sprintf("  Data: n=%d users, p=%d artists, d=%d\n", n, p, d))

# Compute tau
P_temp <- mean(Aselect[upper.tri(Aselect)])
tau <- min(P_temp / 10, 0.001)
cat(sprintf("  Network density: %.4f\n", P_temp))
cat(sprintf("  tau: %.6f\n\n", tau))

# ============================================================================
# RUN ASE
# ============================================================================
cat("Running ASE for initialization...\n")
A_aug <- Aselect
diag(A_aug) <- rowSums(Aselect) / (n - 1)
spec <- eigen(A_aug, symmetric = TRUE)
X_ase <- spec$vectors[, 1:d, drop = FALSE] %*% diag(sqrt(abs(spec$values[1:d])), nrow = d)

# For d=2 with standard RDPG (p=2, q=0)
S <- diag(c(1, 1))
Y_ase <- X_ase %*% S

cat("  ASE complete\n\n")

# ============================================================================
# ANALYZE EDGE PROBABILITIES
# ============================================================================
cat("============================================================================\n")
cat("  EDGE PROBABILITY ANALYSIS\n")
cat("============================================================================\n\n")

cat("Computing all pairwise probabilities P_ij = X_i^T Y_j...\n")

# Compute all edge probabilities
P_matrix <- X_ase %*% t(Y_ase)

# Get upper triangle (unique pairs)
upper_tri_idx <- upper.tri(P_matrix)
P_values <- P_matrix[upper_tri_idx]

cat(sprintf("  Total unique pairs: %d\n\n", length(P_values)))

# Basic statistics
cat("EDGE PROBABILITY STATISTICS:\n")
cat(sprintf("  Min:      %.6f\n", min(P_values)))
cat(sprintf("  Q1:       %.6f\n", quantile(P_values, 0.25)))
cat(sprintf("  Median:   %.6f\n", median(P_values)))
cat(sprintf("  Mean:     %.6f\n", mean(P_values)))
cat(sprintf("  Q3:       %.6f\n", quantile(P_values, 0.75)))
cat(sprintf("  Max:      %.6f\n\n", max(P_values)))

# Check boundary issues
boundary_threshold <- 0.01  # Within 1% of boundary

cat("BOUNDARY ANALYSIS:\n")
cat(sprintf("  Valid range: [%.6f, %.6f]\n", tau, 1-tau))
cat(sprintf("  Boundary threshold: %.4f (1%% of range)\n\n", boundary_threshold))

# Count values near boundaries
below_tau <- sum(P_values < tau)
near_zero <- sum(P_values < boundary_threshold)
above_1_minus_tau <- sum(P_values > 1 - tau)
near_one <- sum(P_values > 1 - boundary_threshold)
in_danger_zone <- sum(P_values < tau + boundary_threshold | P_values > 1 - tau - boundary_threshold)

cat(sprintf("  Below tau (%.6f):         %d (%.2f%%)\n",
            tau, below_tau, 100*below_tau/length(P_values)))
cat(sprintf("  Near 0 (<%.4f):          %d (%.2f%%)\n",
            boundary_threshold, near_zero, 100*near_zero/length(P_values)))
cat(sprintf("  Above 1-tau (%.6f):      %d (%.2f%%)\n",
            1-tau, above_1_minus_tau, 100*above_1_minus_tau/length(P_values)))
cat(sprintf("  Near 1 (>%.4f):          %d (%.2f%%)\n",
            1-boundary_threshold, near_one, 100*near_one/length(P_values)))
cat(sprintf("  In danger zone:           %d (%.2f%%)\n\n",
            in_danger_zone, 100*in_danger_zone/length(P_values)))

# Fisher information quality
cat("FISHER INFORMATION ANALYSIS:\n")
cat("Computing approximate Fisher information...\n")

# For psi(s, tau) = (s-tau)/(1-2*tau), derivative is dpsi/ds = 1/(1-2*tau)
# Fisher weight w = dpsi/ds
# For probabilities near boundaries, Fisher info is poor

# Approximate Fisher weights
fisher_weights <- ifelse(P_values >= tau & P_values <= 1-tau,
                         1 / ((P_values - tau) * (1 - tau - P_values)),
                         Inf)

# Filter out infinite weights
finite_weights <- fisher_weights[is.finite(fisher_weights)]

cat(sprintf("  Pairs with finite weights: %d (%.2f%%)\n",
            length(finite_weights), 100*length(finite_weights)/length(P_values)))

if (length(finite_weights) > 0) {
  cat(sprintf("  Weight statistics (finite only):\n"))
  cat(sprintf("    Min:    %.2e\n", min(finite_weights)))
  cat(sprintf("    Median: %.2e\n", median(finite_weights)))
  cat(sprintf("    Max:    %.2e\n", max(finite_weights)))
  cat(sprintf("    Mean:   %.2e\n\n", mean(finite_weights)))

  # Count poorly-conditioned pairs
  high_weight_threshold <- 1e6
  poorly_conditioned <- sum(finite_weights > high_weight_threshold)
  cat(sprintf("  Poorly conditioned pairs (weight > 1e6): %d (%.2f%%)\n\n",
              poorly_conditioned, 100*poorly_conditioned/length(finite_weights)))
}

# ============================================================================
# FEASIBILITY ASSESSMENT
# ============================================================================
cat("============================================================================\n")
cat("  FEASIBILITY ASSESSMENT\n")
cat("============================================================================\n\n")

# Determine feasibility based on multiple criteria
criteria_passed <- 0
total_criteria <- 5

cat("Checking Fisher-scoring feasibility criteria:\n\n")

# Criterion 1: Most probabilities in valid range
valid_range_pct <- 100 * sum(P_values >= tau & P_values <= 1-tau) / length(P_values)
cat(sprintf("1. Probabilities in [tau, 1-tau]: %.1f%%\n", valid_range_pct))
if (valid_range_pct > 90) {
  cat("   ✓ PASS (>90%%)\n")
  criteria_passed <- criteria_passed + 1
} else {
  cat("   ✗ FAIL (<90%%)\n")
}

# Criterion 2: Few near boundaries
danger_pct <- 100 * in_danger_zone / length(P_values)
cat(sprintf("\n2. Pairs in danger zone: %.1f%%\n", danger_pct))
if (danger_pct < 20) {
  cat("   ✓ PASS (<20%%)\n")
  criteria_passed <- criteria_passed + 1
} else {
  cat("   ✗ FAIL (>20%%)\n")
}

# Criterion 3: Reasonable weight distribution
if (length(finite_weights) > 0) {
  weight_range <- log10(max(finite_weights)) - log10(min(finite_weights))
  cat(sprintf("\n3. Fisher weight range (log10): %.1f\n", weight_range))
  if (weight_range < 10) {
    cat("   ✓ PASS (<10 orders of magnitude)\n")
    criteria_passed <- criteria_passed + 1
  } else {
    cat("   ✗ FAIL (>10 orders of magnitude)\n")
  }
} else {
  cat("\n3. Fisher weights: NO FINITE WEIGHTS\n")
  cat("   ✗ FAIL\n")
}

# Criterion 4: Network density
cat(sprintf("\n4. Network density: %.4f\n", P_temp))
if (P_temp > 0.01) {
  cat("   ✓ PASS (>1%%)\n")
  criteria_passed <- criteria_passed + 1
} else {
  cat("   ⚠ MARGINAL (<1%% - very sparse)\n")
  if (P_temp > 0.001) {
    criteria_passed <- criteria_passed + 0.5
  }
}

# Criterion 5: Condition number estimate
if (length(finite_weights) > 0 && max(finite_weights) < 1e10) {
  cat(sprintf("\n5. Max Fisher weight: %.2e\n", max(finite_weights)))
  if (max(finite_weights) < 1e8) {
    cat("   ✓ PASS (<1e8)\n")
    criteria_passed <- criteria_passed + 1
  } else {
    cat("   ⚠ MARGINAL (1e8-1e10)\n")
    criteria_passed <- criteria_passed + 0.5
  }
} else {
  cat("\n5. Max Fisher weight: EXCESSIVE\n")
  cat("   ✗ FAIL\n")
}

cat("\n============================================================================\n")
cat(sprintf("  CRITERIA PASSED: %.1f / %d\n", criteria_passed, total_criteria))
cat("============================================================================\n\n")

# Final recommendation
if (criteria_passed >= 4) {
  cat("✓ RECOMMENDATION: Fisher-scoring is LIKELY FEASIBLE\n")
  cat("  Expected convergence: Slow but possible\n")
  cat("  Estimated time: 10-30 hours for maxit=30\n\n")
  recommendation <- "FEASIBLE"
} else if (criteria_passed >= 2.5) {
  cat("⚠ RECOMMENDATION: Fisher-scoring is MARGINAL\n")
  cat("  Expected convergence: Very slow or may not converge\n")
  cat("  Estimated time: >24 hours or timeout\n")
  cat("  SUGGESTION: Consider ASE/OSE only, or reduce covariates\n\n")
  recommendation <- "MARGINAL"
} else {
  cat("✗ RECOMMENDATION: Fisher-scoring is NOT FEASIBLE\n")
  cat("  Expected convergence: Will not converge or timeout\n")
  cat("  STRONG SUGGESTION: Use ASE/OSE only for your paper\n\n")
  recommendation <- "NOT_FEASIBLE"
}

# Save diagnostic results
diagnostic_results <- list(
  data_info = list(n = n, p = p, d = d, density = P_temp, tau = tau),
  probability_stats = list(
    min = min(P_values),
    median = median(P_values),
    max = max(P_values),
    below_tau_pct = 100*below_tau/length(P_values),
    above_1minus_tau_pct = 100*above_1_minus_tau/length(P_values),
    danger_zone_pct = danger_pct
  ),
  fisher_info = if(length(finite_weights) > 0) list(
    finite_weight_pct = 100*length(finite_weights)/length(P_values),
    weight_median = median(finite_weights),
    weight_max = max(finite_weights),
    poorly_conditioned_pct = 100*poorly_conditioned/length(finite_weights)
  ) else list(finite_weight_pct = 0),
  criteria_passed = criteria_passed,
  recommendation = recommendation
)

saveRDS(diagnostic_results, "fisher_diagnostic_results.rds")
cat("Diagnostic results saved to: fisher_diagnostic_results.rds\n")

cat("============================================================================\n")
cat("  DIAGNOSTIC COMPLETE\n")
cat("============================================================================\n")
