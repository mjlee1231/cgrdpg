#!/usr/bin/env Rscript
# Design a 3D curve for GRDPG with signature matrix diag(1,1,-1)
# Goal: Wide edge probability range within [0.25, 0.75] for good coverage
# Avoid: Too narrow range (poor coverage) or too close to 0/1 (instability)

library(cgrdpg)

# Test different 3D curve designs
test_curve <- function(n, r, h, c1, c2, c3, name) {
  cat(sprintf("\n=== Testing %s ===\n", name))
  cat(sprintf("Parameters: r=%.3f, h=%.3f, center=(%.3f, %.3f, %.3f)\n", r, h, c1, c2, c3))

  # Generate curve
  theta <- pi * (0:(n-1)) / (n - 1)
  X <- matrix(0, n, 3)
  X[, 1] <- r * cos(theta) + c1
  X[, 2] <- r * sin(theta) + c2
  X[, 3] <- h * theta / pi + c3

  # Compute edge probabilities with signature matrix diag(1,1,-1)
  S <- diag(c(1, 1, -1))
  Y <- X %*% S
  P <- X %*% t(Y)

  # Statistics
  min_p <- min(P)
  max_p <- max(P)
  range_p <- max_p - min_p
  mean_p <- mean(P)

  cat(sprintf("Edge prob range: [%.6f, %.6f] (width=%.4f)\n", min_p, max_p, range_p))
  cat(sprintf("Mean edge prob: %.6f\n", mean_p))

  # Check if all probabilities are in [0.25, 0.75]
  target_min <- 0.25
  target_max <- 0.75
  in_range <- all(P >= target_min & P <= target_max)
  cat(sprintf("All in [%.2f, %.2f]: %s\n", target_min, target_max, in_range))

  if (in_range) {
    cat("✓ THIS DESIGN WORKS!\n")

    # Compute tau
    tau <- min(min_p, 1 - max_p) - 0.001
    cat(sprintf("Suggested tau: %.4f\n", tau))

    # Quality metric: range width (wider is better for coverage)
    cat(sprintf("Quality metric (range width): %.4f\n", range_p))

    return(list(success = TRUE, X = X, P = P, min_p = min_p, max_p = max_p,
                range = range_p, tau = tau))
  } else {
    n_below <- sum(P < target_min)
    n_above <- sum(P > target_max)
    cat(sprintf("✗ Failed: %d below %.2f, %d above %.2f\n", n_below, target_min, n_above, target_max))
    return(list(success = FALSE))
  }
}

n <- 500

cat("============================================================================\n")
cat("  3D CURVE DESIGN for GRDPG with S = diag(1,1,-1)\n")
cat("  Target: WIDE edge probability range within [0.25, 0.75]\n")
cat("  Goal: Balance between coverage quality and numerical stability\n")
cat("============================================================================\n")

# Strategy: Use higher center values to lift the minimum, then increase r and h for spread
# The key insight: c1 and c2 need to be higher to compensate for larger r

# Design 1: Fine-tuned for width ~0.30
result1 <- test_curve(n, r=0.12, h=0.028, c1=0.50, c2=0.50, c3=0.31,
                      name="Width target 0.30 (r=0.12)")

# Design 2: Slightly wider
result2 <- test_curve(n, r=0.13, h=0.030, c1=0.505, c2=0.505, c3=0.31,
                      name="Width target 0.32 (r=0.13)")

# Design 3: Push for 0.35 width
result3 <- test_curve(n, r=0.14, h=0.032, c1=0.51, c2=0.51, c3=0.31,
                      name="Width target 0.35 (r=0.14)")

# Design 4: Maximum safe width attempt
result4 <- test_curve(n, r=0.15, h=0.034, c1=0.515, c2=0.515, c3=0.31,
                      name="Width target 0.37 (r=0.15)")

# Design 5: Alternative with larger h
result5 <- test_curve(n, r=0.13, h=0.035, c1=0.51, c2=0.51, c3=0.315,
                      name="Larger h (r=0.13, h=0.035)")

# Design 6: Conservative but wide
result6 <- test_curve(n, r=0.125, h=0.029, c1=0.505, c2=0.505, c3=0.305,
                      name="Conservative wide (r=0.125)")

# Design 7: Asymmetric centers for variation
result7 <- test_curve(n, r=0.13, h=0.031, c1=0.51, c2=0.505, c3=0.31,
                      name="Asymmetric (r=0.13)")

# Design 8: Push harder
result8 <- test_curve(n, r=0.145, h=0.033, c1=0.515, c2=0.515, c3=0.305,
                      name="Push harder (r=0.145)")

# Design 9: Balance between spread and safety
result9 <- test_curve(n, r=0.135, h=0.031, c1=0.51, c2=0.51, c3=0.308,
                      name="Balanced (r=0.135)")

# Design 10: Maximum attempt with very high centers
result10 <- test_curve(n, r=0.16, h=0.036, c1=0.52, c2=0.52, c3=0.31,
                       name="Max with high centers (r=0.16)")

cat("\n============================================================================\n")
cat("SUMMARY\n")
cat("============================================================================\n")

successful_designs <- list()
if (result1$success) successful_designs[["Design 1"]] <- result1
if (result2$success) successful_designs[["Design 2"]] <- result2
if (result3$success) successful_designs[["Design 3"]] <- result3
if (result4$success) successful_designs[["Design 4"]] <- result4
if (result5$success) successful_designs[["Design 5"]] <- result5
if (result6$success) successful_designs[["Design 6"]] <- result6
if (result7$success) successful_designs[["Design 7"]] <- result7
if (result8$success) successful_designs[["Design 8"]] <- result8
if (result9$success) successful_designs[["Design 9"]] <- result9
if (result10$success) successful_designs[["Design 10"]] <- result10

if (length(successful_designs) > 0) {
  cat(sprintf("Found %d successful design(s):\n\n", length(successful_designs)))

  # Sort by range width (descending)
  ranges <- sapply(successful_designs, function(x) x$range)
  sorted_idx <- order(ranges, decreasing = TRUE)
  sorted_names <- names(successful_designs)[sorted_idx]

  cat("Ranked by edge probability range width (wider = better coverage):\n")
  for (i in seq_along(sorted_names)) {
    name <- sorted_names[i]
    res <- successful_designs[[name]]
    cat(sprintf("  %d. %s: [%.4f, %.4f] (width=%.4f), tau=%.4f\n",
                i, name, res$min_p, res$max_p, res$range, res$tau))
  }

  # Save the best design (widest range)
  best <- successful_designs[[sorted_names[1]]]
  saveRDS(best, "3d_curve_design.rds")
  cat(sprintf("\n*** BEST DESIGN saved to: 3d_curve_design.rds ***\n"))
  cat(sprintf("  Edge prob range: [%.6f, %.6f] (width=%.4f)\n", best$min_p, best$max_p, best$range))
  cat(sprintf("  Recommended tau: %.4f\n", best$tau))
  cat(sprintf("  This design balances:\n"))
  cat(sprintf("    - Wide range for good coverage\n"))
  cat(sprintf("    - Safe distance from 0 and 1 for stability\n"))
} else {
  cat("No successful designs found. Need to adjust parameters.\n")
}
