# Aggregate Diagnostic Oracle Results
# Combines results from diagnostic_2d_tight_oracle and diagnostic_2d_spread_oracle

cat("============================================================================\n")
cat("  Aggregating Diagnostic Oracle Results\n")
cat("============================================================================\n\n")

# Function to aggregate results from individual rep files
aggregate_scenario <- function(pattern, scenario_name) {
  cat(sprintf("\n--- %s ---\n", scenario_name))

  # Find all result files
  files <- list.files(pattern = pattern)

  if (length(files) == 0) {
    cat(sprintf("No files found matching pattern: %s\n", pattern))
    return(NULL)
  }

  cat(sprintf("Found %d result files:\n", length(files)))
  print(files)

  # Read all results
  results_list <- lapply(files, function(f) {
    cat(sprintf("  Reading: %s\n", f))
    readRDS(f)
  })

  # Convert to data frame
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      rep = x$rep,
      SSE = x$SSE,
      fit_time = x$fit_time,
      coverage_true = x$coverage_true,
      coverage_plugin = x$coverage_plugin,
      stringsAsFactors = FALSE
    )
  }))

  # Sort by replication number
  results_df <- results_df[order(results_df$rep), ]

  cat("\nIndividual Results:\n")
  print(results_df, row.names = FALSE)

  # Compute summary statistics
  cat("\n========== SUMMARY STATISTICS ==========\n")
  cat(sprintf("Number of replications: %d\n\n", nrow(results_df)))

  cat(sprintf("SSE:\n"))
  cat(sprintf("  Mean = %.4f\n", mean(results_df$SSE)))
  cat(sprintf("  SD   = %.4f\n", sd(results_df$SSE)))
  cat(sprintf("  Min  = %.4f\n", min(results_df$SSE)))
  cat(sprintf("  Max  = %.4f\n\n", max(results_df$SSE)))

  cat(sprintf("Fit Time (seconds):\n"))
  cat(sprintf("  Mean = %.1f\n", mean(results_df$fit_time)))
  cat(sprintf("  SD   = %.1f\n", sd(results_df$fit_time)))
  cat(sprintf("  Min  = %.1f\n", min(results_df$fit_time)))
  cat(sprintf("  Max  = %.1f\n\n", max(results_df$fit_time)))

  cat(sprintf("Coverage (TRUE G_in):\n"))
  cat(sprintf("  Mean = %.1f%%\n", 100*mean(results_df$coverage_true)))
  cat(sprintf("  SD   = %.1f%%\n", 100*sd(results_df$coverage_true)))
  cat(sprintf("  Min  = %.1f%%\n", 100*min(results_df$coverage_true)))
  cat(sprintf("  Max  = %.1f%%\n\n", 100*max(results_df$coverage_true)))

  cat(sprintf("Coverage (PLUGIN G_in):\n"))
  cat(sprintf("  Mean = %.1f%%\n", 100*mean(results_df$coverage_plugin)))
  cat(sprintf("  SD   = %.1f%%\n", 100*sd(results_df$coverage_plugin)))
  cat(sprintf("  Min  = %.1f%%\n", 100*min(results_df$coverage_plugin)))
  cat(sprintf("  Max  = %.1f%%\n\n", 100*max(results_df$coverage_plugin)))

  cat("========================================\n")

  return(results_df)
}

# Aggregate TIGHT scenario
tight_results <- aggregate_scenario(
  "diagnostic_2d_tight_oracle_rep[0-9]+\\.rds",
  "2D TIGHT Oracle (Adaptive Step Sizes)"
)

# Aggregate SPREAD scenario
spread_results <- aggregate_scenario(
  "diagnostic_2d_spread_oracle_rep[0-9]+\\.rds",
  "2D SPREAD Oracle (Adaptive Step Sizes)"
)

# Save aggregated results
if (!is.null(tight_results)) {
  saveRDS(tight_results, "diagnostic_2d_tight_oracle_aggregated.rds")
  cat("\nTIGHT results saved to: diagnostic_2d_tight_oracle_aggregated.rds\n")
}

if (!is.null(spread_results)) {
  saveRDS(spread_results, "diagnostic_2d_spread_oracle_aggregated.rds")
  cat("SPREAD results saved to: diagnostic_2d_spread_oracle_aggregated.rds\n")
}

# Compare scenarios
if (!is.null(tight_results) && !is.null(spread_results)) {
  cat("\n============================================================================\n")
  cat("  COMPARISON: TIGHT vs SPREAD\n")
  cat("============================================================================\n\n")

  comparison <- data.frame(
    Metric = c("SSE (mean)", "SSE (sd)",
               "Coverage TRUE (mean)", "Coverage TRUE (sd)",
               "Coverage PLUGIN (mean)", "Coverage PLUGIN (sd)",
               "Fit Time (mean)", "Fit Time (sd)"),
    TIGHT = c(
      mean(tight_results$SSE), sd(tight_results$SSE),
      100*mean(tight_results$coverage_true), 100*sd(tight_results$coverage_true),
      100*mean(tight_results$coverage_plugin), 100*sd(tight_results$coverage_plugin),
      mean(tight_results$fit_time), sd(tight_results$fit_time)
    ),
    SPREAD = c(
      mean(spread_results$SSE), sd(spread_results$SSE),
      100*mean(spread_results$coverage_true), 100*sd(spread_results$coverage_true),
      100*mean(spread_results$coverage_plugin), 100*sd(spread_results$coverage_plugin),
      mean(spread_results$fit_time), sd(spread_results$fit_time)
    ),
    stringsAsFactors = FALSE
  )

  print(comparison, row.names = FALSE)
}

cat("\n============================================================================\n")
cat("Aggregation complete!\n")
cat("============================================================================\n")
