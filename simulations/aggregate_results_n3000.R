# Aggregate results from n=3000 array job
library(cgrdpg)

cat("============================================================================\n")
cat("  AGGREGATING RESULTS: n=3000, p_cov=1500\n")
cat("============================================================================\n\n")

# Read all result files
results_dir <- "results_n3000"
result_files <- list.files(results_dir, pattern = "^results_n3000_rep[0-9]+\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  cat("No result files found in", results_dir, "\n")
  quit(status = 1)
}

cat(sprintf("Found %d result files\n\n", length(result_files)))

# Initialize storage
all_results <- list()
coverage_vec <- numeric(length(result_files))
sse_vec <- numeric(length(result_files))
rep_time_vec <- numeric(length(result_files))
fit_time_vec <- numeric(length(result_files))
cov_time_vec <- numeric(length(result_files))
n_inside_vec <- numeric(length(result_files))
n_outside_vec <- numeric(length(result_files))
n_na_vec <- numeric(length(result_files))

# Read all results
for (i in seq_along(result_files)) {
  res <- readRDS(result_files[i])
  all_results[[i]] <- res

  coverage_vec[i] <- res$coverage
  sse_vec[i] <- res$SSE
  rep_time_vec[i] <- res$rep_time
  fit_time_vec[i] <- res$fit_time
  cov_time_vec[i] <- res$cov_time
  n_inside_vec[i] <- res$n_inside
  n_outside_vec[i] <- res$n_outside
  n_na_vec[i] <- res$n_na
}

# Sort by rep_id
rep_ids <- sapply(all_results, function(x) x$rep_id)
sort_idx <- order(rep_ids)
all_results <- all_results[sort_idx]
coverage_vec <- coverage_vec[sort_idx]
sse_vec <- sse_vec[sort_idx]
rep_time_vec <- rep_time_vec[sort_idx]
fit_time_vec <- fit_time_vec[sort_idx]
cov_time_vec <- cov_time_vec[sort_idx]
n_inside_vec <- n_inside_vec[sort_idx]
n_outside_vec <- n_outside_vec[sort_idx]
n_na_vec <- n_na_vec[sort_idx]

# Compute summary statistics
cat("--- COVERAGE STATISTICS ---\n")
cat(sprintf("Mean Coverage:     %.2f%% (SD = %.2f%%)\n",
            100*mean(coverage_vec), 100*sd(coverage_vec)))
cat(sprintf("Median Coverage:   %.2f%%\n", 100*median(coverage_vec)))
cat(sprintf("Min Coverage:      %.2f%% (Rep %d)\n",
            100*min(coverage_vec), which.min(coverage_vec)))
cat(sprintf("Max Coverage:      %.2f%% (Rep %d)\n",
            100*max(coverage_vec), which.max(coverage_vec)))
cat(sprintf("IQR Coverage:      [%.2f%%, %.2f%%]\n",
            100*quantile(coverage_vec, 0.25), 100*quantile(coverage_vec, 0.75)))

cat("\n--- SSE STATISTICS ---\n")
cat(sprintf("Mean SSE:          %.4f (SD = %.4f)\n", mean(sse_vec), sd(sse_vec)))
cat(sprintf("Median SSE:        %.4f\n", median(sse_vec)))
cat(sprintf("Min SSE:           %.4f (Rep %d)\n", min(sse_vec), which.min(sse_vec)))
cat(sprintf("Max SSE:           %.4f (Rep %d)\n", max(sse_vec), which.max(sse_vec)))

cat("\n--- NODES IN/OUT/NA ---\n")
cat(sprintf("Mean Inside:       %.1f / 3000 (SD = %.1f)\n",
            mean(n_inside_vec), sd(n_inside_vec)))
cat(sprintf("Mean Outside:      %.1f / 3000 (SD = %.1f)\n",
            mean(n_outside_vec), sd(n_outside_vec)))
cat(sprintf("Mean NA:           %.1f / 3000 (SD = %.1f)\n",
            mean(n_na_vec), sd(n_na_vec)))

cat("\n--- TIMING STATISTICS ---\n")
cat(sprintf("Mean Total Time:   %.2f hours (SD = %.2f hours)\n",
            mean(rep_time_vec)/3600, sd(rep_time_vec)/3600))
cat(sprintf("Median Total Time: %.2f hours\n", median(rep_time_vec)/3600))
cat(sprintf("Mean Fit Time:     %.2f hours (SD = %.2f hours)\n",
            mean(fit_time_vec)/3600, sd(fit_time_vec)/3600))
cat(sprintf("Mean Cov Time:     %.2f hours (SD = %.2f hours)\n",
            mean(cov_time_vec)/3600, sd(cov_time_vec)/3600))
cat(sprintf("Min Total Time:    %.2f hours (Rep %d)\n",
            min(rep_time_vec)/3600, which.min(rep_time_vec)))
cat(sprintf("Max Total Time:    %.2f hours (Rep %d)\n",
            max(rep_time_vec)/3600, which.max(rep_time_vec)))

# Save aggregated summary
summary_results <- list(
  n_reps = length(result_files),
  coverage_mean = mean(coverage_vec),
  coverage_sd = sd(coverage_vec),
  coverage_median = median(coverage_vec),
  coverage_min = min(coverage_vec),
  coverage_max = max(coverage_vec),
  coverage_q25 = quantile(coverage_vec, 0.25),
  coverage_q75 = quantile(coverage_vec, 0.75),
  sse_mean = mean(sse_vec),
  sse_sd = sd(sse_vec),
  sse_median = median(sse_vec),
  rep_time_mean_hours = mean(rep_time_vec)/3600,
  rep_time_sd_hours = sd(rep_time_vec)/3600,
  fit_time_mean_hours = mean(fit_time_vec)/3600,
  cov_time_mean_hours = mean(cov_time_vec)/3600,
  n = 3000,
  p_cov = 1500,
  d = 3,
  all_results = all_results
)

saveRDS(summary_results, file.path(results_dir, "summary_n3000.rds"))

cat("\n============================================================================\n")
cat(sprintf("Summary saved to: %s\n", file.path(results_dir, "summary_n3000.rds")))
cat("============================================================================\n")
