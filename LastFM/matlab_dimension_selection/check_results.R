#!/usr/bin/env Rscript
# Extract and display MATLAB dimension selection results

library(R.matlab)

cat("============================================================================\n")
cat("MATLAB Dimension Selection Results (maxit=50, tol=0.1)\n")
cat("============================================================================\n\n")

results_dir <- "results_matlab_maxit50_tol0.1"

# Collect all results
all_results <- list()
for (d in 1:5) {
  file_path <- sprintf("%s/results_d%d.mat", results_dir, d)

  if (file.exists(file_path)) {
    mat_data <- readMat(file_path)
    result <- mat_data$result

    all_results[[d]] <- list(
      d = as.integer(result$d),
      neg_log_lik = as.numeric(result$neg.log.lik),
      n_params = as.integer(result$n.params),
      BIC = as.numeric(result$BIC),
      AIC = as.numeric(result$AIC),
      converged = as.logical(result$converged),
      iters = as.integer(result$iters),
      time_sec = as.numeric(result$time.sec)
    )
  }
}

# Display summary table
cat("Summary Table:\n")
cat(sprintf("%-5s | %-15s | %-12s | %-15s | %-15s | %-10s | %-8s | %-12s\n",
            "d", "Neg-Log-Lik", "# Params", "BIC", "AIC", "Converged", "Iters", "Time (hrs)"))
cat(paste(rep("-", 120), collapse=""), "\n")

for (d in 1:5) {
  if (!is.null(all_results[[d]])) {
    r <- all_results[[d]]
    cat(sprintf("%-5d | %15.2f | %12d | %15.2f | %15.2f | %-10s | %8d | %12.2f\n",
                r$d, r$neg_log_lik, r$n_params, r$BIC, r$AIC,
                r$converged, r$iters, r$time_sec/3600))
  }
}

cat("\n")

# Find optimal dimension
BICs <- sapply(all_results, function(r) r$BIC)
AICs <- sapply(all_results, function(r) r$AIC)
optimal_d_BIC <- which.min(BICs)
optimal_d_AIC <- which.min(AICs)

cat("============================================================================\n")
cat("Model Selection:\n")
cat("============================================================================\n")
cat(sprintf("  Optimal dimension (BIC): d = %d (BIC = %.2f)\n", optimal_d_BIC, BICs[optimal_d_BIC]))
cat(sprintf("  Optimal dimension (AIC): d = %d (AIC = %.2f)\n", optimal_d_AIC, AICs[optimal_d_AIC]))
cat("\n")

# Convergence summary
converged_dims <- sapply(all_results, function(r) r$converged)
cat("Convergence Summary:\n")
cat(sprintf("  Converged: %s\n", paste(which(converged_dims), collapse=", ")))
if (any(!converged_dims)) {
  cat(sprintf("  NOT converged: %s\n", paste(which(!converged_dims), collapse=", ")))
}
cat("\n")

# Timing summary
total_time <- sum(sapply(all_results, function(r) r$time_sec)) / 3600
cat("Timing Summary:\n")
for (d in 1:5) {
  if (!is.null(all_results[[d]])) {
    r <- all_results[[d]]
    cat(sprintf("  d=%d: %.2f hours (%d iterations)\n", r$d, r$time_sec/3600, r$iters))
  }
}
cat(sprintf("  Total: %.2f hours\n", total_time))
cat("\n")

cat("============================================================================\n")
