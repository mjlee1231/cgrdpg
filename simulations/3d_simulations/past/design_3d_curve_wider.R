#!/usr/bin/env Rscript
# Design a 3D curve for GRDPG with signature matrix diag(1,1,-1)
# Goal: MAXIMUM edge probability range within [0.25, 0.75] for best coverage

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
                range = range_p, tau = tau, r = r, h = h, c1 = c1, c2 = c2, c3 = c3))
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
cat("  Target: MAXIMUM edge probability range within [0.25, 0.75]\n")
cat("  Strategy: Push to get width as close to 0.50 as possible\n")
cat("============================================================================\n")

# Strategy: Increase r and h even more, with fine-tuned centers
# Current best: r=0.16, h=0.036, width=0.42
# Try to push to width ~0.48-0.49

# Design 1: Increase r slightly
result1 <- test_curve(n, r=0.17, h=0.038, c1=0.525, c2=0.525, c3=0.31,
                      name="r=0.17, h=0.038")

# Design 2: Increase h more
result2 <- test_curve(n, r=0.165, h=0.040, c1=0.525, c2=0.525, c3=0.315,
                      name="r=0.165, h=0.040")

# Design 3: Push r to 0.18
result3 <- test_curve(n, r=0.18, h=0.038, c1=0.53, c2=0.53, c3=0.31,
                      name="r=0.18, h=0.038")

# Design 4: Balanced increase
result4 <- test_curve(n, r=0.175, h=0.040, c1=0.53, c2=0.53, c3=0.315,
                      name="r=0.175, h=0.040")

# Design 5: Maximum r attempt
result5 <- test_curve(n, r=0.19, h=0.038, c1=0.535, c2=0.535, c3=0.31,
                      name="r=0.19, h=0.038")

# Design 6: Maximum h attempt
result6 <- test_curve(n, r=0.17, h=0.042, c1=0.53, c2=0.53, c3=0.32,
                      name="r=0.17, h=0.042")

# Design 7: Both large
result7 <- test_curve(n, r=0.18, h=0.042, c1=0.535, c2=0.535, c3=0.32,
                      name="r=0.18, h=0.042")

# Design 8: Push harder
result8 <- test_curve(n, r=0.185, h=0.041, c1=0.535, c2=0.535, c3=0.315,
                      name="r=0.185, h=0.041")

# Design 9: Conservative max
result9 <- test_curve(n, r=0.175, h=0.041, c1=0.532, c2=0.532, c3=0.316,
                      name="r=0.175, h=0.041")

# Design 10: Fine-tuned
result10 <- test_curve(n, r=0.172, h=0.040, c1=0.528, c2=0.528, c3=0.314,
                       name="r=0.172, h=0.040")

# Design 11: Alternative approach - asymmetric
result11 <- test_curve(n, r=0.18, h=0.040, c1=0.532, c2=0.530, c3=0.315,
                       name="Asymmetric: r=0.18, h=0.040")

# Design 12: Maximum safe attempt
result12 <- test_curve(n, r=0.19, h=0.040, c1=0.537, c2=0.537, c3=0.315,
                       name="Max safe: r=0.19, h=0.040")

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
if (result11$success) successful_designs[["Design 11"]] <- result11
if (result12$success) successful_designs[["Design 12"]] <- result12

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
  cat(sprintf("  Curve parameters:\n"))
  cat(sprintf("    r  = %.3f\n", best$r))
  cat(sprintf("    h  = %.3f\n", best$h))
  cat(sprintf("    c1 = %.3f\n", best$c1))
  cat(sprintf("    c2 = %.3f\n", best$c2))
  cat(sprintf("    c3 = %.3f\n", best$c3))
  cat(sprintf("  Edge prob range: [%.6f, %.6f] (width=%.4f)\n", best$min_p, best$max_p, best$range))
  cat(sprintf("  Recommended tau: %.4f\n", best$tau))
  cat(sprintf("\n  This design maximizes:\n"))
  cat(sprintf("    - Edge probability range for good coverage\n"))
  cat(sprintf("    - While maintaining safe distance from 0 and 1\n"))
} else {
  cat("No successful designs found. Need to adjust parameters.\n")
}
