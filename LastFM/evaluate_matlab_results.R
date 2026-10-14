#!/usr/bin/env Rscript
# Evaluate MATLAB MSLE results using clustering metrics
# Compares MATLAB MSLE with R package results

library(R.matlab)
library(cluster)

# Manual ARI implementation
adjustedRandIndex <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  n <- length(x)
  tab <- table(x, y)
  a <- rowSums(tab)
  b <- colSums(tab)
  sum_comb_tab <- sum(choose(tab, 2))
  sum_comb_a <- sum(choose(a, 2))
  sum_comb_b <- sum(choose(b, 2))
  expected_index <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_index <- (sum_comb_a + sum_comb_b) / 2
  ari <- (sum_comb_tab - expected_index) / (max_index - expected_index)
  return(ari)
}

evaluate_clustering <- function(X, labels, method_name) {
  cat(sprintf("  %s: ", method_name))
  flush.console()

  sil <- silhouette(as.numeric(as.factor(labels)), dist(X))
  avg_sil <- mean(sil[, 3])

  k <- length(unique(labels))
  set.seed(123)
  km <- kmeans(X, centers = k, nstart = 20)
  ari <- adjustedRandIndex(km$cluster, labels)
  bss_tss <- km$betweenss / km$totss

  cat(sprintf("Sil=%.3f, ARI=%.3f, BSS/TSS=%.3f\n", avg_sil, ari, bss_tss))
  flush.console()

  return(list(silhouette = avg_sil, ari = ari, bss_tss = bss_tss))
}

cat("============================================================================\n")
cat("Evaluating MATLAB MSLE Results (maxit=50, tol=0.1)\n")
cat("============================================================================\n\n")

# Load data for labels
cat("Loading data for labels...\n")
load("CovariatesMatrix.Rdata")
X <- as.matrix(sparseX)
edgelist <- read.csv("lastfm_asia_edges.csv") + 1
label <- read.csv('lastfm_asia_target.csv')$target + 1

# Same preprocessing as in dimension_selection_single.m
ind5 <- which(label != 5)
X <- X[ind5,]
label <- label[ind5]
X <- X[, colSums(X) > 0]

labelnew <- label
label[labelnew > 5] <- label[labelnew > 5] - 1

sizes <- summary(as.factor(label))
class_select <- names(sizes)[sizes > 300 & sizes < 1000]
ind.select <- which(label %in% class_select)
Xselect <- X[ind.select,]
labelselect <- label[ind.select]

dartist_all <- colSums(X)
dartist <- colSums(Xselect)
prop <- dartist / dartist_all
n_select <- nrow(Xselect)
prob <- 1 - round(min(n_select/2, 600)/ncol(Xselect), 4)
Xselect <- Xselect[, prop > quantile(prop, probs=prob, na.rm=TRUE)]

dnode <- rowSums(Xselect)
like <- which(dnode > 0)
Xselect <- Xselect[like,]
labelselect <- labelselect[like]

# Remove low-degree nodes (matching MATLAB)
A <- matrix(0, nrow(X), nrow(X))
for (i in 1:nrow(edgelist)) A[edgelist[i,1], edgelist[i,2]] <- 1
A <- A + t(A)
A <- A[ind5, ind5]
Aselect <- A[ind.select, ind.select]
Aselect <- Aselect[like, like]

dselect <- rowSums(Aselect)
while(sum(dselect <= 1) > 0) {
  keep <- which(dselect > 1)
  Aselect <- Aselect[keep,keep]
  Xselect <- Xselect[keep,]
  labelselect <- labelselect[keep]
  dselect <- rowSums(Aselect)
}

cat(sprintf("Final data: n=%d, k=%d countries\n\n", nrow(Xselect), length(unique(labelselect))))

# Process each dimension
for (d in c(2, 3)) {
  mat_file <- sprintf("matlab_dimension_selection/results_matlab_maxit50_tol0.1/results_d%d.mat", d)

  if (!file.exists(mat_file)) {
    cat(sprintf("Skipping d=%d: file not found\n\n", d))
    next
  }

  cat("============================================================================\n")
  cat(sprintf("d = %d\n", d))
  cat("============================================================================\n\n")

  # Load MATLAB results
  mat_data <- readMat(mat_file)
  result <- mat_data$result

  # Extract latent positions
  X_matlab <- result[[10]]  # latent.positions is 10th element

  cat(sprintf("MATLAB MSLE latent positions: %d x %d\n", nrow(X_matlab), ncol(X_matlab)))
  cat(sprintf("  Converged: %s\n", ifelse(result[[6]][1,1] == 1, "YES", "NO")))
  cat(sprintf("  Iterations: %d\n", as.integer(result[[7]][1,1])))
  cat(sprintf("  Time: %.2f hours\n", result[[8]][1,1]/3600))
  cat(sprintf("  BIC: %.2f\n\n", result[[4]][1,1]))

  # Evaluate clustering
  cat("Clustering Performance:\n")
  perf_matlab <- evaluate_clustering(X_matlab, labelselect, "MATLAB MSLE")

  # Compare with R results if available
  r_result_file <- sprintf("d%d/results/lastfm_d%d_complete_results.rds", d, d)

  if (file.exists(r_result_file)) {
    cat("\n")
    r_result <- readRDS(r_result_file)

    cat("R Package Results (for comparison):\n")
    cat(sprintf("  R MSLE:  Sil=%.3f, ARI=%.3f, BSS/TSS=%.3f (%.2fh, %s)\n",
                r_result$performance$fisher$silhouette,
                r_result$performance$fisher$ari,
                r_result$performance$fisher$bss_tss,
                r_result$timing$fisher/3600,
                ifelse(r_result$fisher_info$converged, "converged", "NOT converged")))

    cat(sprintf("  R ASE:   Sil=%.3f, ARI=%.3f, BSS/TSS=%.3f\n",
                r_result$performance$ase$silhouette,
                r_result$performance$ase$ari,
                r_result$performance$ase$bss_tss))

    cat(sprintf("  R OSE:   Sil=%.3f, ARI=%.3f, BSS/TSS=%.3f\n",
                r_result$performance$ose$silhouette,
                r_result$performance$ose$ari,
                r_result$performance$ose$bss_tss))

    cat("\nComparison (MATLAB vs R MSLE):\n")
    cat(sprintf("  Silhouette: %.3f vs %.3f (diff: %+.3f)\n",
                perf_matlab$silhouette,
                r_result$performance$fisher$silhouette,
                perf_matlab$silhouette - r_result$performance$fisher$silhouette))

    cat(sprintf("  ARI:        %.3f vs %.3f (diff: %+.3f)\n",
                perf_matlab$ari,
                r_result$performance$fisher$ari,
                perf_matlab$ari - r_result$performance$fisher$ari))

    cat(sprintf("  BSS/TSS:    %.3f vs %.3f (diff: %+.3f)\n",
                perf_matlab$bss_tss,
                r_result$performance$fisher$bss_tss,
                perf_matlab$bss_tss - r_result$performance$fisher$bss_tss))

    cat(sprintf("  Time:       %.2fh vs %.2fh (speedup: %.1fx)\n",
                result[[8]][1,1]/3600,
                r_result$timing$fisher/3600,
                r_result$timing$fisher / result[[8]][1,1]))
  } else {
    cat(sprintf("\n(R results not available for d=%d)\n", d))
  }

  cat("\n")
}

cat("============================================================================\n")
cat("Evaluation Complete\n")
cat("============================================================================\n")
