# Aggregate results from all replications (n=800)
# Run this after all array jobs complete

cat("============================================================================\n")
cat("  Aggregating Results from n=800 Simulation\n")
cat("============================================================================\n\n")

# Find all result files in the results_n800 directory
result_dir <- "results_n800"
if (!dir.exists(result_dir)) {
  stop("Directory 'results_n800' not found! Make sure replications have completed.")
}

result_files <- list.files(path = result_dir,
                           pattern = "^results_n800_rep[0-9]+\\.rds$",
                           full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in results_n800/! Make sure replications have completed.")
}

cat(sprintf("Found %d result files:\n", length(result_files)))
for (f in result_files) {
  cat(sprintf("  - %s\n", f))
}
cat("\n")

# Read all results
all_results <- lapply(result_files, readRDS)

# Extract statistics
rep_ids <- sapply(all_results, function(x) x$rep_id)
coverages <- sapply(all_results, function(x) x$coverage)
SSEs <- sapply(all_results, function(x) x$SSE)
rep_times <- sapply(all_results, function(x) x$rep_time)
fit_times <- sapply(all_results, function(x) x$fit_time)
n_nas <- sapply(all_results, function(x) x$n_na)

n <- all_results[[1]]$n
p_cov <- all_results[[1]]$p_cov

cat("============================================================================\n")
cat("  AGGREGATE RESULTS\n")
cat("============================================================================\n\n")

cat(sprintf("Scenario: n=%d, p_cov=%d\n", n, p_cov))
cat(sprintf("Number of replications: %d\n\n", length(all_results)))

cat("--- Per-Replication Results ---\n")
cat(sprintf("%-5s %-12s %-10s %-10s %-8s\n", "Rep", "Coverage", "SSE", "Time(min)", "NAs"))
cat(sprintf("%-5s %-12s %-10s %-10s %-8s\n", "---", "--------", "------", "---------", "---"))
for (i in 1:length(all_results)) {
  cat(sprintf("%-5d %-12.1f%% %-10.4f %-10.1f %-8d\n",
              rep_ids[i], 100*coverages[i], SSEs[i], rep_times[i]/60, n_nas[i]))
}

cat("\n--- Coverage Statistics ---\n")
cat(sprintf("Mean coverage:   %.1f%% (Expected: 95.0%%)\n", 100*mean(coverages)))
cat(sprintf("SD coverage:     %.1f%%\n", 100*sd(coverages)))
cat(sprintf("Min coverage:    %.1f%%\n", 100*min(coverages)))
cat(sprintf("Max coverage:    %.1f%%\n", 100*max(coverages)))
cat(sprintf("Median coverage: %.1f%%\n\n", 100*median(coverages)))

cat("--- SSE Statistics ---\n")
cat(sprintf("Mean SSE:        %.4f\n", mean(SSEs)))
cat(sprintf("SD SSE:          %.4f\n", sd(SSEs)))
cat(sprintf("Min SSE:         %.4f\n", min(SSEs)))
cat(sprintf("Max SSE:         %.4f\n", max(SSEs)))
cat(sprintf("Median SSE:      %.4f\n\n", median(SSEs)))

cat("--- Timing Statistics ---\n")
cat(sprintf("Mean time per rep:    %.1f minutes (%.2f hours)\n",
            mean(rep_times)/60, mean(rep_times)/3600))
cat(sprintf("Total wall time:      %.1f minutes (%.2f hours)\n",
            max(rep_times)/60, max(rep_times)/3600))
cat(sprintf("  (All ran in parallel, limited by slowest)\n"))
cat(sprintf("Mean fitting time:    %.1f minutes\n", mean(fit_times)/60))
cat(sprintf("Fitting time range:   %.1f - %.1f minutes\n\n",
            min(fit_times)/60, max(fit_times)/60))

cat("--- Numerical Issues ---\n")
cat(sprintf("Total NAs across all reps:  %d\n", sum(n_nas)))
cat(sprintf("Mean NAs per rep:           %.1f\n", mean(n_nas)))
cat(sprintf("Max NAs in a single rep:    %d\n\n", max(n_nas)))

# Save aggregate summary
summary <- list(
  n = n,
  p_cov = p_cov,
  n_reps = length(all_results),
  mean_coverage = mean(coverages),
  sd_coverage = sd(coverages),
  mean_SSE = mean(SSEs),
  sd_SSE = sd(SSEs),
  mean_time = mean(rep_times),
  coverages = coverages,
  SSEs = SSEs,
  rep_times = rep_times
)

summary_file <- file.path(result_dir, "aggregate_summary_n800.rds")
saveRDS(summary, summary_file)
cat(sprintf("Aggregate summary saved to: %s\n", summary_file))

cat("============================================================================\n")
