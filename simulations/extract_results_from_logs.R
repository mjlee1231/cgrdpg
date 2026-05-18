# Extract Results from Log Files
# When .rds files failed to save, extract results from .out log files

extract_from_log <- function(log_file) {
  if (!file.exists(log_file)) {
    cat(sprintf("Log file not found: %s\n", log_file))
    return(NULL)
  }

  cat(sprintf("Reading: %s\n", log_file))
  lines <- readLines(log_file)

  # Find the results table line
  result_line <- grep("^[0-9]+\\s+[0-9]+\\.", lines, value = TRUE)

  if (length(result_line) == 0) {
    cat("  No results found in log\n")
    return(NULL)
  }

  # Parse the line: "1     2.9732     1290.6     89.0            90.4"
  parts <- strsplit(trimws(result_line[1]), "\\s+")[[1]]

  result <- list(
    rep = as.integer(parts[1]),
    SSE = as.numeric(parts[2]),
    fit_time = as.numeric(parts[3]),
    coverage_true = as.numeric(parts[4]) / 100,
    coverage_plugin = as.numeric(parts[5]) / 100
  )

  cat(sprintf("  Rep %d: SSE=%.4f, Time=%.1fs, CovTRUE=%.1f%%, CovPLUGIN=%.1f%%\n",
              result$rep, result$SSE, result$fit_time,
              100*result$coverage_true, 100*result$coverage_plugin))

  return(result)
}

# Extract SPREAD results
cat("============================================================================\n")
cat("  Extracting SPREAD Oracle Results from Log Files\n")
cat("============================================================================\n\n")

spread_logs <- list.files("logs", pattern = "diag_spread_rep.*\\.out", full.names = TRUE)
cat(sprintf("Found %d log files\n\n", length(spread_logs)))

spread_results <- lapply(spread_logs, extract_from_log)
spread_results <- spread_results[!sapply(spread_results, is.null)]

if (length(spread_results) > 0) {
  spread_df <- do.call(rbind, lapply(spread_results, function(x) {
    data.frame(
      rep = x$rep,
      SSE = x$SSE,
      fit_time = x$fit_time,
      coverage_true = x$coverage_true,
      coverage_plugin = x$coverage_plugin,
      stringsAsFactors = FALSE
    )
  }))

  spread_df <- spread_df[order(spread_df$rep), ]

  cat("\n============================================================================\n")
  cat("  SPREAD ORACLE RESULTS (ALL REPLICATIONS)\n")
  cat("============================================================================\n\n")

  print(spread_df, row.names = FALSE)

  cat("\n========== SUMMARY STATISTICS ==========\n")
  cat(sprintf("Number of replications: %d\n\n", nrow(spread_df)))

  cat(sprintf("SSE:\n"))
  cat(sprintf("  Mean = %.4f\n", mean(spread_df$SSE)))
  cat(sprintf("  SD   = %.4f\n", sd(spread_df$SSE)))
  cat(sprintf("  Min  = %.4f\n", min(spread_df$SSE)))
  cat(sprintf("  Max  = %.4f\n\n", max(spread_df$SSE)))

  cat(sprintf("Fit Time (seconds):\n"))
  cat(sprintf("  Mean = %.1f\n", mean(spread_df$fit_time)))
  cat(sprintf("  SD   = %.1f\n", sd(spread_df$fit_time)))
  cat(sprintf("  Min  = %.1f\n", min(spread_df$fit_time)))
  cat(sprintf("  Max  = %.1f\n\n", max(spread_df$fit_time)))

  cat(sprintf("Coverage (TRUE G_in):\n"))
  cat(sprintf("  Mean = %.1f%%\n", 100*mean(spread_df$coverage_true)))
  cat(sprintf("  SD   = %.1f%%\n", 100*sd(spread_df$coverage_true)))
  cat(sprintf("  Min  = %.1f%%\n", 100*min(spread_df$coverage_true)))
  cat(sprintf("  Max  = %.1f%%\n\n", 100*max(spread_df$coverage_true)))

  cat(sprintf("Coverage (PLUGIN G_in):\n"))
  cat(sprintf("  Mean = %.1f%%\n", 100*mean(spread_df$coverage_plugin)))
  cat(sprintf("  SD   = %.1f%%\n", 100*sd(spread_df$coverage_plugin)))
  cat(sprintf("  Min  = %.1f%%\n", 100*min(spread_df$coverage_plugin)))
  cat(sprintf("  Max  = %.1f%%\n\n", 100*max(spread_df$coverage_plugin)))

  cat("========================================\n\n")

  # Save extracted results
  saveRDS(spread_df, "diagnostic_2d_spread_oracle_aggregated.rds")
  cat("Results saved to: diagnostic_2d_spread_oracle_aggregated.rds\n")

} else {
  cat("\nNo results could be extracted from log files.\n")
}

cat("\n============================================================================\n")
