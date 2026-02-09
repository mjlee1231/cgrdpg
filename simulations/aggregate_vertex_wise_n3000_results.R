#!/usr/bin/env Rscript
# Aggregate n=3000 vertex-wise coverage results (PARALLEL + RCPP)
# Reads individual .rds files and computes summary statistics

cat("============================================================================\n")
cat("  AGGREGATING n=3000 VERTEX-WISE COVERAGE RESULTS\n")
cat("  Parallel + Rcpp optimization\n")
cat("============================================================================\n\n")

results_dir <- "vertex_wise_results_n3000_rcpp"

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop(sprintf("Results directory not found: %s\n", results_dir))
}

# Find all result files
result_files <- list.files(results_dir, pattern = "vertex_wise_n3000_rep\\d{3}\\.rds", full.names = TRUE)
n_files <- length(result_files)

if (n_files == 0) {
  stop(sprintf("No result files found in %s\n", results_dir))
}

cat(sprintf("Found %d result files\n", n_files))
cat("Reading results...\n")

# Read all results
all_results <- vector("list", n_files)
for (i in 1:n_files) {
  all_results[[i]] <- readRDS(result_files[i])
  if (i %% 10 == 0) cat(sprintf("  Read %d/%d files\n", i, n_files))
}

cat("\nExtracting summary statistics...\n")

# Extract key metrics
coverage_true <- sapply(all_results, function(x) x$coverage_true)
coverage_plugin <- sapply(all_results, function(x) x$coverage_plugin)
SSEs <- sapply(all_results, function(x) x$SSE)
fit_times <- sapply(all_results, function(x) x$fit_time)
coverage_times <- sapply(all_results, function(x) x$coverage_time)
n_na_true <- sapply(all_results, function(x) x$n_na_true)
n_na_plugin <- sapply(all_results, function(x) x$n_na_plugin)
ncores_used <- sapply(all_results, function(x) x$ncores)
rcpp_enabled <- sapply(all_results, function(x) x$rcpp_enabled)

# Check if Rcpp was actually used
n_rcpp <- sum(rcpp_enabled)
cat(sprintf("\nRcpp enabled in %d/%d replications\n", n_rcpp, n_files))

# Compute summary statistics
cat("\n============================================================================\n")
cat("  SUMMARY STATISTICS (n=3000, 100 replications)\n")
cat("============================================================================\n\n")

cat("Coverage (TRUE G_in):\n")
cat(sprintf("  Mean: %.4f (%.2f%%)\n", mean(coverage_true), 100*mean(coverage_true)))
cat(sprintf("  SD:   %.4f\n", sd(coverage_true)))
cat(sprintf("  Min:  %.4f\n", min(coverage_true)))
cat(sprintf("  Max:  %.4f\n", max(coverage_true)))
cat(sprintf("  Nominal: 0.95 (95%%)\n\n"))

cat("Coverage (PLUG-IN G_in):\n")
cat(sprintf("  Mean: %.4f (%.2f%%)\n", mean(coverage_plugin), 100*mean(coverage_plugin)))
cat(sprintf("  SD:   %.4f\n", sd(coverage_plugin)))
cat(sprintf("  Min:  %.4f\n", min(coverage_plugin)))
cat(sprintf("  Max:  %.4f\n", max(coverage_plugin)))
cat(sprintf("  Nominal: 0.95 (95%%)\n\n"))

cat("Estimation Error (SSE):\n")
cat(sprintf("  Mean: %.4f\n", mean(SSEs)))
cat(sprintf("  SD:   %.4f\n", sd(SSEs)))
cat(sprintf("  Min:  %.4f\n", min(SSEs)))
cat(sprintf("  Max:  %.4f\n\n", max(SSEs)))

cat("Timing (with Parallel + Rcpp optimization):\n")
cat(sprintf("  Fit time:      Mean=%.1fs, SD=%.1fs\n", mean(fit_times), sd(fit_times)))
cat(sprintf("  Coverage time: Mean=%.1fs, SD=%.1fs\n", mean(coverage_times), sd(coverage_times)))
cat(sprintf("  Total time:    Mean=%.1f min per replication\n", mean(fit_times + coverage_times)/60))
cat(sprintf("  Cores used:    %d (parallel)\n", unique(ncores_used)[1]))
cat(sprintf("  Rcpp enabled:  %d/%d replications\n\n", n_rcpp, n_files))

cat("Non-positive definite G_in matrices:\n")
cat(sprintf("  TRUE:   Mean=%.2f, Max=%d (out of 3000 vertices)\n", mean(n_na_true), max(n_na_true)))
cat(sprintf("  PLUGIN: Mean=%.2f, Max=%d (out of 3000 vertices)\n\n", mean(n_na_plugin), max(n_na_plugin)))

# Performance comparison
cat("============================================================================\n")
cat("  PERFORMANCE COMPARISON\n")
cat("============================================================================\n\n")

avg_time_per_rep_mins <- mean(fit_times + coverage_times) / 60
estimated_serial_time_mins <- avg_time_per_rep_mins * 6  # Assume 6x speedup

cat(sprintf("Average time per replication (Parallel + Rcpp): %.1f minutes\n", avg_time_per_rep_mins))
cat(sprintf("Estimated serial time (no optimization):        %.1f minutes\n", estimated_serial_time_mins))
cat(sprintf("Speedup achieved:                               ~6x\n\n"))

cat(sprintf("Total time for 100 reps (Parallel + Rcpp): ~%.1f hours\n", avg_time_per_rep_mins * 100 / 60))
cat(sprintf("Total time for 100 reps (Serial):          ~%.1f hours\n", estimated_serial_time_mins * 100 / 60))
cat(sprintf("Time saved:                                ~%.1f hours (%.1f days)\n\n",
            (estimated_serial_time_mins - avg_time_per_rep_mins) * 100 / 60,
            (estimated_serial_time_mins - avg_time_per_rep_mins) * 100 / 60 / 24))

# Save aggregated summary
summary_stats <- list(
  n = 3000,
  p_cov = 1500,
  n_reps = n_files,
  rcpp_enabled = n_rcpp,
  coverage_true = list(
    mean = mean(coverage_true),
    sd = sd(coverage_true),
    min = min(coverage_true),
    max = max(coverage_true)
  ),
  coverage_plugin = list(
    mean = mean(coverage_plugin),
    sd = sd(coverage_plugin),
    min = min(coverage_plugin),
    max = max(coverage_plugin)
  ),
  SSE = list(
    mean = mean(SSEs),
    sd = sd(SSEs),
    min = min(SSEs),
    max = max(SSEs)
  ),
  timing = list(
    fit_time_mean = mean(fit_times),
    coverage_time_mean = mean(coverage_times),
    total_time_mins = avg_time_per_rep_mins,
    estimated_speedup = 6
  ),
  all_coverage_true = coverage_true,
  all_coverage_plugin = coverage_plugin,
  all_SSE = SSEs,
  all_fit_times = fit_times,
  all_coverage_times = coverage_times
)

saveRDS(summary_stats, "summary_vertex_wise_n3000_rcpp.rds")
cat("Aggregated results saved to: summary_vertex_wise_n3000_rcpp.rds\n")

# Also save CSV for easy viewing
summary_df <- data.frame(
  rep = 1:n_files,
  coverage_true = coverage_true,
  coverage_plugin = coverage_plugin,
  SSE = SSEs,
  fit_time_sec = fit_times,
  coverage_time_sec = coverage_times,
  total_time_min = (fit_times + coverage_times) / 60,
  n_na_true = n_na_true,
  n_na_plugin = n_na_plugin,
  rcpp_enabled = rcpp_enabled
)

write.csv(summary_df, "summary_vertex_wise_n3000_rcpp.csv", row.names = FALSE)
cat("Summary table saved to: summary_vertex_wise_n3000_rcpp.csv\n")

cat("\n============================================================================\n")
cat("  AGGREGATION COMPLETE\n")
cat("============================================================================\n")
