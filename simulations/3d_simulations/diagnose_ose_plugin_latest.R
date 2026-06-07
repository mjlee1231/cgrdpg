#!/usr/bin/env Rscript
# Diagnostic: Check OSE plugin precision matrices from latest run

library(cgrdpg)

# Load the result
res <- readRDS("results_3d_ase_ose_cgrdpg_n1000/rep_001.rds")

cat("=== OSE Plugin Diagnostic ===\n\n")
cat(sprintf("OSE SSE: %.4f\n", res$sse["ose"]))
cat(sprintf("OSE true coverage: %.1f%%\n", 100 * res$overall_cov["ose_true"]))
cat(sprintf("OSE plugin coverage: %.1f%%\n\n", 100 * res$overall_cov["ose_plugin"]))

# Check how many NAs in plugin
cat(sprintf("Number of NA values in ose_plugin: %d\n\n", res$n_na["ose_plugin"]))

cat("This suggests precision matrices are degenerate (eigenvalue < 1e-10)\n")
cat("despite OSE SSE improving to 293.9.\n\n")
cat("OSE estimation is still not good enough for reliable precision computation.\n")
