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

# Function to process scenario results
process_scenario <- function(pattern, scenario_name, output_file) {
  cat("============================================================================\n")
  cat(sprintf("  Extracting %s Oracle Results from Log Files\n", scenario_name))
  cat("============================================================================\n\n")

  log_files <- list.files("logs", pattern = pattern, full.names = TRUE)
  cat(sprintf("Found %d log files\n\n", length(log_files)))

  if (length(log_files) == 0) {
    cat(sprintf("No log files found matching pattern: %s\n\n", pattern))
    return(NULL)
  }

  results <- lapply(log_files, extract_from_log)
  results <- results[!sapply(results, is.null)]

  if (length(results) == 0) {
    cat("\nNo results could be extracted from log files.\n\n")
    return(NULL)
  }

  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      rep = x$rep,
      SSE = x$SSE,
      fit_time = x$fit_time,
      coverage_true = x$coverage_true,
      coverage_plugin = x$coverage_plugin,
      stringsAsFactors = FALSE
    )
  }))

  results_df <- results_df[order(results_df$rep), ]

  cat("\n============================================================================\n")
  cat(sprintf("  %s ORACLE RESULTS (ALL REPLICATIONS)\n", toupper(scenario_name)))
  cat("============================================================================\n\n")

  print(results_df, row.names = FALSE)

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

  cat("========================================\n\n")

  # Save extracted results
  saveRDS(results_df, output_file)
  cat(sprintf("Results saved to: %s\n\n", output_file))

  return(results_df)
}

# Extract TIGHT results
tight_df <- process_scenario(
  "diag_tight_rep.*\\.out",
  "TIGHT",
  "diagnostic_2d_tight_oracle_aggregated.rds"
)

# Extract SPREAD results
spread_df <- process_scenario(
  "diag_spread_rep.*\\.out",
  "SPREAD",
  "diagnostic_2d_spread_oracle_aggregated.rds"
)

# Compare scenarios if both exist
if (!is.null(tight_df) && !is.null(spread_df)) {
  cat("============================================================================\n")
  cat("  COMPARISON: TIGHT vs SPREAD\n")
  cat("============================================================================\n\n")

  comparison <- data.frame(
    Metric = c("SSE (mean)", "SSE (sd)",
               "Coverage TRUE (mean)", "Coverage TRUE (sd)",
               "Coverage PLUGIN (mean)", "Coverage PLUGIN (sd)",
               "Fit Time (mean)", "Fit Time (sd)"),
    TIGHT = c(
      mean(tight_df$SSE), sd(tight_df$SSE),
      100*mean(tight_df$coverage_true), 100*sd(tight_df$coverage_true),
      100*mean(tight_df$coverage_plugin), 100*sd(tight_df$coverage_plugin),
      mean(tight_df$fit_time), sd(tight_df$fit_time)
    ),
    SPREAD = c(
      mean(spread_df$SSE), sd(spread_df$SSE),
      100*mean(spread_df$coverage_true), 100*sd(spread_df$coverage_true),
      100*mean(spread_df$coverage_plugin), 100*sd(spread_df$coverage_plugin),
      mean(spread_df$fit_time), sd(spread_df$fit_time)
    ),
    stringsAsFactors = FALSE
  )

  print(comparison, row.names = FALSE)
  cat("\n")
}

cat("============================================================================\n")
cat("Extraction complete!\n")
cat("============================================================================\n")
