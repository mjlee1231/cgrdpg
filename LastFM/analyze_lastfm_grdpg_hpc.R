#!/usr/bin/env Rscript
# LastFM Analysis with GRDPG: Compare Fisher vs ASE vs OSE
# OPTIMIZED FOR HPC: maxit=3 for faster execution
# Research questions:
# 1. Fit GRDPG model with covariates
# 2. Compare Fisher-scoring vs ASE vs OSE
# 3. Compute vertex-wise confidence regions
# 4. Evaluate alignment with country labels (community detection)

library(cgrdpg)
library(igraph)
library(Matrix)
library(irlba)
library(ggplot2)
library(cluster)

# Manual implementation of Adjusted Rand Index
adjustedRandIndex <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  n <- length(x)

  # Contingency table
  tab <- table(x, y)

  # Row and column sums
  a <- rowSums(tab)
  b <- colSums(tab)

  # Compute index
  sum_comb_tab <- sum(choose(tab, 2))
  sum_comb_a <- sum(choose(a, 2))
  sum_comb_b <- sum(choose(b, 2))
  expected_index <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_index <- (sum_comb_a + sum_comb_b) / 2

  ari <- (sum_comb_tab - expected_index) / (max_index - expected_index)
  return(ari)
}

cat("============================================================================\n")
cat("  LASTFM ASIA SOCIAL NETWORK ANALYSIS WITH GRDPG (HPC VERSION)\n")
cat("  Comparing: FISHER (with covariates) vs ASE vs OSE (no covariates)\n")
cat("============================================================================\n\n")

# ============================================================================
# STEP 1: DATA PREPROCESSING (following dimension_reduction.R)
# ============================================================================
cat("STEP 1: Loading and preprocessing data...\n")

# Load covariate matrix
load("CovariatesMatrix.Rdata")
X <- as.matrix(sparseX)
n <- dim(X)[1]
p <- dim(X)[2]
cat(sprintf("  Initial dimensions: n=%d users, p=%d artists\n", n, p))

# Load adjacency matrix
edgelist <- read.csv("lastfm_asia_edges.csv")
edgelist <- edgelist + 1  # Convert 0-indexed to 1-indexed
A <- matrix(0, n, n)
for (i in 1:dim(edgelist)[1]) {
  A[edgelist[i, 1], edgelist[i, 2]] <- 1
}
A <- A + t(A)  # Symmetrize

# Load country labels
label <- read.csv('lastfm_asia_target.csv')
label <- label$target + 1  # Convert 0-indexed to 1-indexed

cat(sprintf("  Countries before filtering: %d\n", length(unique(label))))

# Delete country with only 17 users (country 5)
ind5 <- which(label != 5)
A <- A[ind5, ind5]
X <- X[ind5, ]
label <- label[ind5]
X <- X[, colSums(X) > 0]  # Remove artists with no listeners

# Relabel countries after removal
labelnew <- label
label[labelnew > 5] <- label[labelnew > 5] - 1

n <- dim(A)[1]
p <- dim(X)[2]
dartist_all <- colSums(X)

cat(sprintf("\n  After removing country 5: n=%d, p=%d\n", n, p))

# Filter by country size: 300-1000 users
sizes <- summary(as.factor(label))
class_select <- names(sizes)[sizes > 300 & sizes < 1000]

cat(sprintf("  Selected countries (300-1000 users): %s\n", paste(class_select, collapse=", ")))

# Select corresponding nodes
ind.select <- which(label %in% as.vector(class_select))
Aselect <- A[ind.select, ind.select]
Xselect <- X[ind.select, ]
labelselect <- label[ind.select]

n <- dim(Aselect)[1]
p <- dim(Xselect)[2]
kk <- length(unique(labelselect))

cat(sprintf("  After country selection: n=%d, p=%d, k=%d countries\n", n, p, kk))

# Select regionally popular artists
dartist <- colSums(Xselect)
prop <- dartist / dartist_all
prob <- 1 - round(min(n/2, 600)/p, 4)
Xselect <- Xselect[, prop > quantile(prop, probs = prob, na.rm = TRUE)]

cat(sprintf("  After artist selection: p=%d artists\n", ncol(Xselect)))

# Remove users who like no artists
dnode <- rowSums(Xselect)
like <- which(dnode > 0)
Xselect <- Xselect[like, ]
Aselect <- Aselect[like, like]
labelselect <- labelselect[like]

# Remove low-degree nodes
dselect <- rowSums(Aselect)
tau_degree <- 1
while (sum(dselect <= tau_degree) > 0) {
  keep_idx <- which(dselect > tau_degree)
  Aselect <- Aselect[keep_idx, keep_idx]
  Xselect <- Xselect[keep_idx, ]
  labelselect <- labelselect[keep_idx]
  dselect <- rowSums(Aselect)
}

n <- dim(Aselect)[1]
p <- dim(Xselect)[2]

cat(sprintf("\n  FINAL DATA: n=%d users, p=%d artists, k=%d countries\n", n, p, length(unique(labelselect))))
cat("  Final country distribution:\n")
print(table(labelselect))

# Transpose covariate matrix for cgrdpg (needs p × n)
B <- t(Xselect)
cat(sprintf("  Covariate matrix B: %d × %d (p × n)\n\n", nrow(B), ncol(B)))

# ============================================================================
# STEP 2: FIT MODELS
# ============================================================================
cat("============================================================================\n")
cat("STEP 2: Fitting models...\n")
cat("============================================================================\n\n")

# Determine tau for Fisher
P_temp <- mean(Aselect[upper.tri(Aselect)])
tau <- min(P_temp / 10, 0.001)
cat(sprintf("  Network density: %.4f\n", P_temp))
cat(sprintf("  Selected tau: %.6f\n\n", tau))

# Determine d from eigenvalue gap analysis
# (Largest gap at position 2, suggesting d=2)
d <- 2
cat(sprintf("  Using d=%d (selected from eigenvalue gap analysis)\n\n", d))

# METHOD 1: FISHER-SCORING (with covariates)
cat("--- METHOD 1: FISHER-SCORING (with covariates) ---\n")
fisher_start <- Sys.time()
fit_fisher <- fit_grdpg_cov(Aselect, B, d = d, p = d, q = 0,
                             maxit = 30, tol = 0.01, tau = tau)
fisher_time <- as.numeric(difftime(Sys.time(), fisher_start, units = "secs"))
cat(sprintf("  Fitting time: %.2f seconds (%.2f minutes)\n", fisher_time, fisher_time/60))
cat(sprintf("  Converged: %s\n", ifelse(fit_fisher$converged, "YES", "NO")))
cat(sprintf("  Final loss: %.4f\n\n", tail(fit_fisher$loss, 1)))

X_fisher <- fit_fisher$X

# METHOD 2: ASE (no covariates)
cat("--- METHOD 2: ASE (no covariates) ---\n")
ase_start <- Sys.time()
A_aug <- Aselect
diag(A_aug) <- rowSums(Aselect) / (n - 1)
spec <- eigen(A_aug, symmetric = TRUE)
X_ase <- spec$vectors[, 1:d, drop = FALSE] %*% diag(sqrt(abs(spec$values[1:d])), nrow = d)
ase_time <- as.numeric(difftime(Sys.time(), ase_start, units = "secs"))
cat(sprintf("  Fitting time: %.2f seconds\n\n", ase_time))

# METHOD 3: OSE (no covariates, 1-step from ASE)
cat("--- METHOD 3: OSE (no covariates, 1-step refinement) ---\n")
ose_start <- Sys.time()

# One-step Newton-Raphson update from ASE
X_ose <- X_ase  # Start from ASE
eps_clip <- 0.001

for (i in 1:n) {
  x_i <- X_ase[i, ]
  idx_j <- setdiff(1:n, i)

  # Compute probabilities
  p_ij <- pmax(pmin(X_ase[idx_j, ] %*% x_i, 1 - eps_clip), eps_clip)

  # Residuals
  resid <- Aselect[i, idx_j] - p_ij

  # Weights
  w <- 1 / (p_ij * (1 - p_ij))

  # Gradient and Hessian
  grad <- t(X_ase[idx_j, ]) %*% (resid * w)
  H <- t(X_ase[idx_j, ]) %*% (X_ase[idx_j, ] * as.vector(w))

  # Update
  X_ose[i, ] <- x_i + solve(H + diag(1e-6, d)) %*% grad
}

ose_time <- as.numeric(difftime(Sys.time(), ose_start, units = "secs"))
cat(sprintf("  Fitting time: %.2f seconds\n\n", ose_time))

# ============================================================================
# STEP 3: VERTEX-WISE CONFIDENCE REGIONS
# ============================================================================
cat("============================================================================\n")
cat("STEP 3: Computing vertex-wise confidence regions...\n")
cat("============================================================================\n\n")

# For demonstration, compute confidence region size for a sample of vertices
sample_size <- min(100, n)
sample_idx <- sample(1:n, sample_size)

cat(sprintf("  Computing confidence ellipses for %d sample vertices...\n", sample_size))

# Helper function to compute Fisher information for vertex i
compute_fisher_precision <- function(i, X, tau) {
  n <- nrow(X)
  d <- ncol(X)
  x_i <- X[i, ]

  # Network part
  s <- X %*% x_i
  w <- dpsi(s, tau = tau)

  G_net <- matrix(0, d, d)
  for (j in 1:n) {
    G_net <- G_net + w[j] * outer(X[j, ], X[j, ])
  }

  # Covariate part (simplified - assumes Z contribution)
  G_cov <- diag(ncol(B) / n, d)

  G <- (G_net + G_cov) / (n + ncol(B))
  return(G)
}

# Compute ellipse sizes (determinant of precision matrix)
fisher_ellipse_sizes <- numeric(sample_size)
ase_ellipse_sizes <- numeric(sample_size)
ose_ellipse_sizes <- numeric(sample_size)

for (idx in 1:sample_size) {
  i <- sample_idx[idx]

  # Fisher
  G_fisher <- compute_fisher_precision(i, X_fisher, tau)
  fisher_ellipse_sizes[idx] <- det(G_fisher)

  # ASE (simplified sandwich variance)
  Delta <- t(X_ase) %*% X_ase
  fisher_ellipse_sizes_ase <- det(Delta) / n^2
  ase_ellipse_sizes[idx] <- fisher_ellipse_sizes_ase

  # OSE (Fisher information)
  x_i <- X_ose[i, ]
  idx_j <- setdiff(1:n, i)
  p_ij <- pmax(pmin(X_ose[idx_j, ] %*% x_i, eps_clip), 1 - eps_clip)
  w <- 1 / (p_ij * (1 - p_ij))
  G_ose <- t(X_ose[idx_j, ]) %*% (X_ose[idx_j, ] * as.vector(w)) / n
  ose_ellipse_sizes[idx] <- det(G_ose)
}

cat(sprintf("  Average ellipse size (det of precision):\n"))
cat(sprintf("    FISHER: %.6f\n", mean(fisher_ellipse_sizes, na.rm=TRUE)))
cat(sprintf("    ASE:    %.6f\n", mean(ase_ellipse_sizes, na.rm=TRUE)))
cat(sprintf("    OSE:    %.6f\n\n", mean(ose_ellipse_sizes, na.rm=TRUE)))

# ============================================================================
# STEP 4: COMMUNITY DETECTION PERFORMANCE
# ============================================================================
cat("============================================================================\n")
cat("STEP 4: Evaluating alignment with country labels...\n")
cat("============================================================================\n\n")

# Performance metrics: ARI, Silhouette score, k-means clustering accuracy

evaluate_clustering <- function(X, labels, method_name) {
  cat(sprintf("--- %s ---\n", method_name))

  # Silhouette score
  sil <- silhouette(as.numeric(as.factor(labels)), dist(X))
  avg_sil <- mean(sil[, 3])
  cat(sprintf("  Silhouette score: %.4f\n", avg_sil))

  # K-means clustering (k = number of countries)
  k <- length(unique(labels))
  set.seed(123)
  km <- kmeans(X, centers = k, nstart = 20)

  # Adjusted Rand Index
  ari <- adjustedRandIndex(km$cluster, labels)
  cat(sprintf("  ARI (k-means vs true labels): %.4f\n", ari))

  # Between-cluster sum of squares / Total sum of squares
  bss_tss <- km$betweenss / km$totss
  cat(sprintf("  BSS/TSS ratio: %.4f\n", bss_tss))

  # Within-cluster sum of squares (average)
  avg_withinss <- mean(km$withinss)
  cat(sprintf("  Avg within-cluster SS: %.2f\n\n", avg_withinss))

  return(list(
    silhouette = avg_sil,
    ari = ari,
    bss_tss = bss_tss
  ))
}

perf_fisher <- evaluate_clustering(X_fisher, labelselect, "FISHER (with covariates)")
perf_ase <- evaluate_clustering(X_ase, labelselect, "ASE (no covariates)")
perf_ose <- evaluate_clustering(X_ose, labelselect, "OSE (no covariates)")

# ============================================================================
# SUMMARY
# ============================================================================
cat("============================================================================\n")
cat("  SUMMARY: METHOD COMPARISON\n")
cat("============================================================================\n\n")

cat("TIMING:\n")
cat(sprintf("  FISHER: %.2f sec (%.2f min)\n", fisher_time, fisher_time/60))
cat(sprintf("  ASE:    %.2f sec\n", ase_time))
cat(sprintf("  OSE:    %.2f sec\n\n", ose_time))

cat("COUNTRY LABEL ALIGNMENT (higher is better):\n")
cat(sprintf("  Method          Silhouette    ARI      BSS/TSS\n"))
cat(sprintf("  %-15s %.4f        %.4f   %.4f\n", "FISHER", perf_fisher$silhouette, perf_fisher$ari, perf_fisher$bss_tss))
cat(sprintf("  %-15s %.4f        %.4f   %.4f\n", "ASE", perf_ase$silhouette, perf_ase$ari, perf_ase$bss_tss))
cat(sprintf("  %-15s %.4f        %.4f   %.4f\n\n", "OSE", perf_ose$silhouette, perf_ose$ari, perf_ose$bss_tss))

cat("INTERPRETATION:\n")
cat("  - Silhouette: How well-separated clusters are (-1 to 1, higher better)\n")
cat("  - ARI: Adjusted Rand Index with true labels (0 to 1, higher better)\n")
cat("  - BSS/TSS: Between/Total sum of squares (0 to 1, higher better)\n\n")

# Determine winner
scores <- data.frame(
  Method = c("FISHER", "ASE", "OSE"),
  Silhouette = c(perf_fisher$silhouette, perf_ase$silhouette, perf_ose$silhouette),
  ARI = c(perf_fisher$ari, perf_ase$ari, perf_ose$ari),
  BSS_TSS = c(perf_fisher$bss_tss, perf_ase$bss_tss, perf_ose$bss_tss)
)

cat("RANKING BY METRIC:\n")
for (metric in c("Silhouette", "ARI", "BSS_TSS")) {
  ranked <- scores[order(-scores[[metric]]), ]
  cat(sprintf("  %s: %s\n", metric, paste(ranked$Method, collapse=" > ")))
}

# Save results
results <- list(
  data = list(
    n = n,
    p = p,
    k = length(unique(labelselect)),
    labels = labelselect,
    A = Aselect,
    B = B
  ),
  latent_positions = list(
    fisher = X_fisher,
    ase = X_ase,
    ose = X_ose
  ),
  performance = list(
    fisher = perf_fisher,
    ase = perf_ase,
    ose = perf_ose
  ),
  timing = list(
    fisher = fisher_time,
    ase = ase_time,
    ose = ose_time
  )
)

saveRDS(results, "lastfm_grdpg_analysis_results_hpc.rds")
cat("\n\nResults saved to: lastfm_grdpg_analysis_results_hpc.rds\n")

cat("============================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("============================================================================\n")
