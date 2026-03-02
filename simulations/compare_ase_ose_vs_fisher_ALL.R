#!/usr/bin/env Rscript
# Compare ASE/OSE (no covariates) vs Fisher-scoring (with covariates)
# Across n=500, 1000, 2000, 3000

cat("============================================================================\n")
cat("  COMPREHENSIVE COMPARISON: ASE/OSE vs FISHER-SCORING\n")
cat("  Comparing coverage rates across sample sizes (n=500, 1000, 2000, 3000)\n")
cat("============================================================================\n\n")

# --- Load ASE/OSE results (no covariates) ---
cat("Loading ASE/OSE results (no covariates)...\n")
ase_ose_500 <- readRDS("vertexwise_coverage_ase_ose_n500.rds")
ase_ose_1000 <- readRDS("vertexwise_coverage_ase_ose_n1000.rds")
ase_ose_2000 <- readRDS("vertexwise_coverage_ase_ose_n2000.rds")
ase_ose_3000 <- readRDS("vertexwise_coverage_ase_ose_n3000.rds")

# --- Load Fisher-scoring results (with covariates) ---
cat("Loading Fisher-scoring results (with covariates)...\n")

# n=500 - optimized version
fisher_500 <- readRDS("vertex_wise_results_n500_optimized/vertex_wise_coverage_aggregated.rds")

# n=1000 - try multiple possible locations
if (file.exists("vertex_wise_results_n1000/vertex_wise_coverage_aggregated_n1000.rds")) {
  fisher_1000 <- readRDS("vertex_wise_results_n1000/vertex_wise_coverage_aggregated_n1000.rds")
} else if (file.exists("vertex_wise_results_n1000/vertex_wise_coverage_aggregated.rds")) {
  fisher_1000 <- readRDS("vertex_wise_results_n1000/vertex_wise_coverage_aggregated.rds")
} else {
  stop("Cannot find n=1000 Fisher-scoring aggregated results")
}

# n=2000
if (file.exists("vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated_n2000.rds")) {
  fisher_2000 <- readRDS("vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated_n2000.rds")
} else if (file.exists("vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated.rds")) {
  fisher_2000 <- readRDS("vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated.rds")
} else {
  stop("Cannot find n=2000 Fisher-scoring aggregated results")
}

# n=3000
if (file.exists("summary_vertex_wise_n3000_rcpp.rds")) {
  fisher_3000 <- readRDS("summary_vertex_wise_n3000_rcpp.rds")
} else if (file.exists("vertexwise_coverage_n3000.rds")) {
  fisher_3000 <- readRDS("vertexwise_coverage_n3000.rds")
} else {
  stop("Cannot find n=3000 Fisher-scoring aggregated results")
}

cat("\n")

# --- Helper function to extract mean coverage ---
extract_coverage <- function(result_obj, method = "fisher") {
  if (method == "ase_ose") {
    list(
      ase_true = mean(result_obj$vertex_coverage$ase_true, na.rm = TRUE),
      ase_plugin = mean(result_obj$vertex_coverage$ase_plugin, na.rm = TRUE),
      ose_true = mean(result_obj$vertex_coverage$ose_true, na.rm = TRUE),
      ose_plugin = mean(result_obj$vertex_coverage$ose_plugin, na.rm = TRUE)
    )
  } else {
    # Fisher-scoring - try multiple possible structures
    fisher_true <- NA
    fisher_plugin <- NA

    # Try structure 1: vertex_coverage_mean/vertex_coverage_plugin_mean
    if (!is.null(result_obj$vertex_coverage_mean) && !is.na(result_obj$vertex_coverage_mean)) {
      fisher_true <- result_obj$vertex_coverage_mean
      fisher_plugin <- result_obj$vertex_coverage_plugin_mean
    }
    # Try structure 2: coverage_true/coverage_plugin with $mean
    else if (!is.null(result_obj$coverage_true) && is.list(result_obj$coverage_true)) {
      fisher_true <- result_obj$coverage_true$mean
      fisher_plugin <- result_obj$coverage_plugin$mean
    }
    # Try structure 3: direct vectors all_coverage_true/all_coverage_plugin
    else if (!is.null(result_obj$all_coverage_true) && is.numeric(result_obj$all_coverage_true)) {
      fisher_true <- mean(result_obj$all_coverage_true, na.rm = TRUE)
      fisher_plugin <- mean(result_obj$all_coverage_plugin, na.rm = TRUE)
    }
    # Try structure 4: vertex_coverage with vectors
    else if (!is.null(result_obj$vertex_coverage) && is.list(result_obj$vertex_coverage)) {
      if (!is.null(result_obj$vertex_coverage$coverage_true)) {
        fisher_true <- mean(result_obj$vertex_coverage$coverage_true, na.rm = TRUE)
        fisher_plugin <- mean(result_obj$vertex_coverage$coverage_plugin, na.rm = TRUE)
      }
    }
    # Try structure 5: Direct vectors vertex_coverage_true/vertex_coverage_plugin
    else if (!is.null(result_obj$vertex_coverage_true) && is.numeric(result_obj$vertex_coverage_true)) {
      fisher_true <- mean(result_obj$vertex_coverage_true, na.rm = TRUE)
      fisher_plugin <- mean(result_obj$vertex_coverage_plugin, na.rm = TRUE)
    }

    if (is.na(fisher_true) || is.na(fisher_plugin)) {
      cat("Warning: Could not extract coverage from Fisher result\n")
      cat("Available names:\n")
      print(names(result_obj))
      if (!is.null(result_obj$vertex_coverage)) {
        cat("\nvertex_coverage names:\n")
        print(names(result_obj$vertex_coverage))
      }
    }

    list(
      fisher_true = fisher_true,
      fisher_plugin = fisher_plugin
    )
  }
}

# --- Extract mean coverage for all scenarios ---
cat("Extracting coverage statistics...\n\n")

cov_ase_500 <- extract_coverage(ase_ose_500, "ase_ose")
cov_ase_1000 <- extract_coverage(ase_ose_1000, "ase_ose")
cov_ase_2000 <- extract_coverage(ase_ose_2000, "ase_ose")
cov_ase_3000 <- extract_coverage(ase_ose_3000, "ase_ose")

cov_fisher_500 <- extract_coverage(fisher_500, "fisher")
cov_fisher_1000 <- extract_coverage(fisher_1000, "fisher")
cov_fisher_2000 <- extract_coverage(fisher_2000, "fisher")
cov_fisher_3000 <- extract_coverage(fisher_3000, "fisher")

# --- Create comparison table ---
cat("============================================================================\n")
cat("  VERTEX-WISE COVERAGE COMPARISON (Mean across vertices and replications)\n")
cat("============================================================================\n\n")

cat("LEGEND:\n")
cat("  ASE-TRUE:       ASE with true precision matrix (no covariates)\n")
cat("  ASE-PLUGIN:     ASE with plug-in precision matrix (no covariates)\n")
cat("  OSE-TRUE:       OSE with true precision matrix (no covariates)\n")
cat("  OSE-PLUGIN:     OSE with plug-in precision matrix (no covariates)\n")
cat("  FISHER-TRUE:    Fisher-scoring with true precision matrix (with covariates)\n")
cat("  FISHER-PLUGIN:  Fisher-scoring with plug-in precision (with covariates)\n")
cat("  Nominal:        95% (target coverage)\n\n")

# Print comparison table
comparison_df <- data.frame(
  n = c(500, 1000, 2000, 3000),
  ASE_TRUE = c(cov_ase_500$ase_true, cov_ase_1000$ase_true, cov_ase_2000$ase_true, cov_ase_3000$ase_true),
  ASE_PLUGIN = c(cov_ase_500$ase_plugin, cov_ase_1000$ase_plugin, cov_ase_2000$ase_plugin, cov_ase_3000$ase_plugin),
  OSE_TRUE = c(cov_ase_500$ose_true, cov_ase_1000$ose_true, cov_ase_2000$ose_true, cov_ase_3000$ose_true),
  OSE_PLUGIN = c(cov_ase_500$ose_plugin, cov_ase_1000$ose_plugin, cov_ase_2000$ose_plugin, cov_ase_3000$ose_plugin),
  FISHER_TRUE = c(cov_fisher_500$fisher_true, cov_fisher_1000$fisher_true, cov_fisher_2000$fisher_true, cov_fisher_3000$fisher_true),
  FISHER_PLUGIN = c(cov_fisher_500$fisher_plugin, cov_fisher_1000$fisher_plugin, cov_fisher_2000$fisher_plugin, cov_fisher_3000$fisher_plugin)
)

cat("Coverage Rates (proportion):\n")
print(comparison_df, digits = 4, row.names = FALSE)
cat("\n")

cat("Coverage Rates (percentage):\n")
comparison_pct <- comparison_df
comparison_pct[, -1] <- comparison_pct[, -1] * 100
print(comparison_pct, digits = 3, row.names = FALSE)
cat("\n")

# --- Key comparisons ---
cat("============================================================================\n")
cat("  KEY COMPARISONS\n")
cat("============================================================================\n\n")

cat("1. ASE vs OSE (both with TRUE precision, no covariates):\n")
for (i in 1:4) {
  n_val <- comparison_df$n[i]
  ase <- comparison_df$ASE_TRUE[i]
  ose <- comparison_df$OSE_TRUE[i]
  diff <- (ose - ase) * 100
  cat(sprintf("   n=%d: OSE=%.2f%%, ASE=%.2f%%, Difference=%.2f%%\n",
              n_val, ose*100, ase*100, diff))
}
cat("\n")

cat("2. Fisher-scoring vs OSE (both with TRUE precision):\n")
cat("   (Fisher has covariates, OSE does not)\n")
for (i in 1:4) {
  n_val <- comparison_df$n[i]
  fisher <- comparison_df$FISHER_TRUE[i]
  ose <- comparison_df$OSE_TRUE[i]
  diff <- (fisher - ose) * 100
  cat(sprintf("   n=%d: Fisher=%.2f%%, OSE=%.2f%%, Difference=%.2f%%\n",
              n_val, fisher*100, ose*100, diff))
}
cat("\n")

cat("3. TRUE vs PLUGIN precision (for OSE):\n")
for (i in 1:4) {
  n_val <- comparison_df$n[i]
  true_p <- comparison_df$OSE_TRUE[i]
  plugin <- comparison_df$OSE_PLUGIN[i]
  diff <- (true_p - plugin) * 100
  cat(sprintf("   n=%d: TRUE=%.2f%%, PLUGIN=%.2f%%, Difference=%.2f%%\n",
              n_val, true_p*100, plugin*100, diff))
}
cat("\n")

cat("4. TRUE vs PLUGIN precision (for Fisher-scoring):\n")
for (i in 1:4) {
  n_val <- comparison_df$n[i]
  true_p <- comparison_df$FISHER_TRUE[i]
  plugin <- comparison_df$FISHER_PLUGIN[i]
  diff <- (true_p - plugin) * 100
  cat(sprintf("   n=%d: TRUE=%.2f%%, PLUGIN=%.2f%%, Difference=%.2f%%\n",
              n_val, true_p*100, plugin*100, diff))
}
cat("\n")

# --- Coverage convergence to nominal ---
cat("============================================================================\n")
cat("  CONVERGENCE TO NOMINAL 95% COVERAGE\n")
cat("============================================================================\n\n")

cat("Distance from nominal (95%) - TRUE precision methods:\n")
for (i in 1:4) {
  n_val <- comparison_df$n[i]
  ase_dist <- abs(comparison_df$ASE_TRUE[i] - 0.95) * 100
  ose_dist <- abs(comparison_df$OSE_TRUE[i] - 0.95) * 100
  fisher_dist <- abs(comparison_df$FISHER_TRUE[i] - 0.95) * 100
  cat(sprintf("   n=%d: ASE=%.2f%%, OSE=%.2f%%, Fisher=%.2f%%\n",
              n_val, ase_dist, ose_dist, fisher_dist))
}
cat("\n")

# --- SSE Comparison (if available) ---
cat("============================================================================\n")
cat("  SSE COMPARISON (ASE vs OSE, no covariates)\n")
cat("============================================================================\n\n")

extract_sse <- function(result_obj) {
  list(
    ase_mean = mean(result_obj$sse$ase),
    ose_mean = mean(result_obj$sse$ose),
    improvement = 100 * (mean(result_obj$sse$ase) - mean(result_obj$sse$ose)) / mean(result_obj$sse$ase)
  )
}

sse_500 <- extract_sse(ase_ose_500)
sse_1000 <- extract_sse(ase_ose_1000)
sse_2000 <- extract_sse(ase_ose_2000)
sse_3000 <- extract_sse(ase_ose_3000)

cat("Mean SSE across replications:\n")
cat(sprintf("n=500:  ASE=%.4f, OSE=%.4f, Improvement=%.2f%%\n",
            sse_500$ase_mean, sse_500$ose_mean, sse_500$improvement))
cat(sprintf("n=1000: ASE=%.4f, OSE=%.4f, Improvement=%.2f%%\n",
            sse_1000$ase_mean, sse_1000$ose_mean, sse_1000$improvement))
cat(sprintf("n=2000: ASE=%.4f, OSE=%.4f, Improvement=%.2f%%\n",
            sse_2000$ase_mean, sse_2000$ose_mean, sse_2000$improvement))
cat(sprintf("n=3000: ASE=%.4f, OSE=%.4f, Improvement=%.2f%%\n\n",
            sse_3000$ase_mean, sse_3000$ose_mean, sse_3000$improvement))

# --- Save comparison results ---
comparison_results <- list(
  comparison_table = comparison_df,
  sse_comparison = data.frame(
    n = c(500, 1000, 2000, 3000),
    ase_sse = c(sse_500$ase_mean, sse_1000$ase_mean, sse_2000$ase_mean, sse_3000$ase_mean),
    ose_sse = c(sse_500$ose_mean, sse_1000$ose_mean, sse_2000$ose_mean, sse_3000$ose_mean),
    ose_improvement_pct = c(sse_500$improvement, sse_1000$improvement, sse_2000$improvement, sse_3000$improvement)
  ),
  raw_data = list(
    ase_ose_500 = ase_ose_500,
    ase_ose_1000 = ase_ose_1000,
    ase_ose_2000 = ase_ose_2000,
    ase_ose_3000 = ase_ose_3000,
    fisher_500 = fisher_500,
    fisher_1000 = fisher_1000,
    fisher_2000 = fisher_2000,
    fisher_3000 = fisher_3000
  )
)

saveRDS(comparison_results, "comparison_ase_ose_vs_fisher_ALL.rds")
write.csv(comparison_df, "comparison_ase_ose_vs_fisher_ALL.csv", row.names = FALSE)

cat("Comparison results saved to:\n")
cat("  - comparison_ase_ose_vs_fisher_ALL.rds\n")
cat("  - comparison_ase_ose_vs_fisher_ALL.csv\n")

cat("\n============================================================================\n")
cat("  COMPARISON COMPLETE\n")
cat("============================================================================\n")
