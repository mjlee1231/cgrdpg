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
    tryCatch({
      mat_data <- readMat(file_path)
      result <- mat_data$result

      # Handle nested structure - extract values safely
      all_results[[d]] <- list(
        d = as.integer(result[[1]]),  # First element is d
        neg_log_lik = as.numeric(result[[2]]),  # neg.log.lik
        n_params = as.integer(result[[3]]),  # n.params
        BIC = as.numeric(result[[4]]),  # BIC
        AIC = as.numeric(result[[5]]),  # AIC
        converged = as.logical(result[[6]]),  # converged
        iters = as.integer(result[[7]]),  # iters
        time_sec = as.numeric(result[[8]])  # time.sec
      )
    }, error = function(e) {
      cat(sprintf("Warning: Could not read %s: %s\n", file_path, e$message))
      all_results[[d]] <<- NULL
    })
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

# Find optimal dimension (filter out NULLs)
valid_results <- all_results[!sapply(all_results, is.null)]
if (length(valid_results) > 0) {
  BICs <- sapply(valid_results, function(r) r$BIC)
  AICs <- sapply(valid_results, function(r) r$AIC)
  dims <- sapply(valid_results, function(r) r$d)

  optimal_idx_BIC <- which.min(BICs)
  optimal_idx_AIC <- which.min(AICs)

  cat("============================================================================\n")
  cat("Model Selection:\n")
  cat("============================================================================\n")
  cat(sprintf("  Optimal dimension (BIC): d = %d (BIC = %.2f)\n",
              dims[optimal_idx_BIC], BICs[optimal_idx_BIC]))
  cat(sprintf("  Optimal dimension (AIC): d = %d (AIC = %.2f)\n",
              dims[optimal_idx_AIC], AICs[optimal_idx_AIC]))
  cat("\n")

  # Convergence summary
  converged_dims <- sapply(valid_results, function(r) as.logical(r$converged))
  converged_d <- dims[converged_dims]
  not_converged_d <- dims[!converged_dims]

  cat("Convergence Summary:\n")
  if (length(converged_d) > 0) {
    cat(sprintf("  Converged: d = %s\n", paste(converged_d, collapse=", ")))
  }
  if (length(not_converged_d) > 0) {
    cat(sprintf("  NOT converged: d = %s\n", paste(not_converged_d, collapse=", ")))
  }
  cat("\n")
} else {
  cat("\nNo valid results found!\n\n")
}

# Timing summary
if (length(valid_results) > 0) {
  total_time <- sum(sapply(valid_results, function(r) r$time_sec)) / 3600
  cat("Timing Summary:\n")
  for (r in valid_results) {
    cat(sprintf("  d=%d: %.2f hours (%d iterations)\n", r$d, r$time_sec/3600, r$iters))
  }
  cat(sprintf("  Total: %.2f hours\n", total_time))
  cat("\n")
}

cat("============================================================================\n")
