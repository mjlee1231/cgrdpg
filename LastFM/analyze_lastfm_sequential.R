#!/usr/bin/env Rscript
# LastFM Analysis: Sequential approach with early results saving
# Strategy: Run ASE/OSE first (fast), save results, then Fisher (slow)
# This ensures we have ASE/OSE results even if Fisher times out

library(cgrdpg)
library(igraph)
library(Matrix)
library(irlba)
library(ggplot2)
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

cat("============================================================================\n")
cat("  LASTFM SEQUENTIAL ANALYSIS: ASE/OSE first, then Fisher\n")
cat("============================================================================\n\n")

# ============================================================================
# PREPROCESSING
# ============================================================================
cat("STEP 1: Loading and preprocessing data...\n")

load("CovariatesMatrix.Rdata")
X <- as.matrix(sparseX)
edgelist <- read.csv("lastfm_asia_edges.csv") + 1
A <- matrix(0, nrow(X), nrow(X))
for (i in 1:nrow(edgelist)) A[edgelist[i,1], edgelist[i,2]] <- 1
A <- A + t(A)

label <- read.csv('lastfm_asia_target.csv')$target + 1
ind5 <- which(label != 5)
A <- A[ind5, ind5]; X <- X[ind5,]; label <- label[ind5]
X <- X[, colSums(X) > 0]

labelnew <- label; label[labelnew > 5] <- label[labelnew > 5] - 1

sizes <- summary(as.factor(label))
class_select <- names(sizes)[sizes > 300 & sizes < 1000]
ind.select <- which(label %in% class_select)
Aselect <- A[ind.select, ind.select]
Xselect <- X[ind.select,]
labelselect <- label[ind.select]

dartist_all <- colSums(X)
dartist <- colSums(Xselect)
prop <- dartist / dartist_all
prob <- 1 - round(min(nrow(Aselect)/2, 600)/ncol(Xselect), 4)
Xselect <- Xselect[, prop > quantile(prop, probs=prob, na.rm=TRUE)]

dnode <- rowSums(Xselect); like <- which(dnode > 0)
Xselect <- Xselect[like,]; Aselect <- Aselect[like,like]; labelselect <- labelselect[like]

dselect <- rowSums(Aselect)
while(sum(dselect <= 1) > 0) {
  keep <- which(dselect > 1)
  Aselect <- Aselect[keep,keep]; Xselect <- Xselect[keep,]; labelselect <- labelselect[keep]
  dselect <- rowSums(Aselect)
}

n <- nrow(Aselect); p <- ncol(Xselect)
B <- t(Xselect)
d <- 2
eps_clip <- 0.001

cat(sprintf("  FINAL: n=%d, p=%d, k=%d countries, d=%d\n\n",
            n, p, length(unique(labelselect)), d))

# ============================================================================
# PART 1: ASE and OSE (FAST - Complete in <1 minute)
# ============================================================================
cat("============================================================================\n")
cat("PART 1: ASE and OSE (fast methods)\n")
cat("============================================================================\n\n")

# ASE
cat("Running ASE...\n")
ase_start <- Sys.time()
A_aug <- Aselect
diag(A_aug) <- rowSums(Aselect) / (n - 1)
spec <- eigen(A_aug, symmetric = TRUE)
X_ase <- spec$vectors[, 1:d, drop = FALSE] %*% diag(sqrt(abs(spec$values[1:d])), nrow = d)
ase_time <- as.numeric(difftime(Sys.time(), ase_start, units = "secs"))
cat(sprintf("  ASE completed in %.2f seconds\n\n", ase_time))

# OSE
cat("Running OSE (from ASE)...\n")
ose_start <- Sys.time()
X_ose <- X_ase

for (i in 1:n) {
  x_i <- X_ase[i, ]
  idx_j <- setdiff(1:n, i)
  p_ij <- pmax(pmin(X_ase[idx_j, ] %*% x_i, 1 - eps_clip), eps_clip)
  resid <- Aselect[i, idx_j] - p_ij
  w <- 1 / (p_ij * (1 - p_ij))
  grad <- t(X_ase[idx_j, ]) %*% (resid * w)
  H <- t(X_ase[idx_j, ]) %*% (X_ase[idx_j, ] * as.vector(w))
  X_ose[i, ] <- x_i + solve(H + diag(1e-6, d)) %*% grad
}
ose_time <- as.numeric(difftime(Sys.time(), ose_start, units = "secs"))
cat(sprintf("  OSE completed in %.2f seconds\n\n", ose_time))

# Evaluate ASE and OSE
cat("Evaluating ASE and OSE...\n")
eval_start <- Sys.time()

evaluate_clustering <- function(X, labels, method_name) {
  cat(sprintf("  %s: ", method_name))
  sil <- silhouette(as.numeric(as.factor(labels)), dist(X))
  avg_sil <- mean(sil[, 3])
  k <- length(unique(labels))
  set.seed(123)
  km <- kmeans(X, centers = k, nstart = 20)
  ari <- adjustedRandIndex(km$cluster, labels)
  bss_tss <- km$betweenss / km$totss
  cat(sprintf("Sil=%.3f, ARI=%.3f, BSS/TSS=%.3f\n", avg_sil, ari, bss_tss))
  return(list(silhouette = avg_sil, ari = ari, bss_tss = bss_tss))
}

perf_ase <- evaluate_clustering(X_ase, labelselect, "ASE")
perf_ose <- evaluate_clustering(X_ose, labelselect, "OSE")

eval_time <- as.numeric(difftime(Sys.time(), eval_start, units = "secs"))

# SAVE ASE/OSE RESULTS NOW (before Fisher)
ase_ose_results <- list(
  data = list(n = n, p = p, k = length(unique(labelselect)), labels = labelselect),
  latent_positions = list(ase = X_ase, ose = X_ose),
  performance = list(ase = perf_ase, ose = perf_ose),
  timing = list(ase = ase_time, ose = ose_time, eval = eval_time)
)

saveRDS(ase_ose_results, "lastfm_ase_ose_results.rds")
cat(sprintf("\n✓ ASE/OSE results saved to: lastfm_ase_ose_results.rds\n"))
cat(sprintf("  Total time so far: %.1f seconds\n\n", ase_time + ose_time + eval_time))

# ============================================================================
# PART 2: FISHER-SCORING (SLOW - May take many hours)
# ============================================================================
cat("============================================================================\n")
cat("PART 2: Fisher-scoring (slow, may take hours)\n")
cat("============================================================================\n\n")

P_temp <- mean(Aselect[upper.tri(Aselect)])
tau <- min(P_temp / 10, 0.001)

cat(sprintf("Starting Fisher-scoring with maxit=30...\n"))
cat(sprintf("Network density: %.4f, tau: %.6f\n", P_temp, tau))
cat(sprintf("This may take 10-20 hours. Progress is buffered.\n\n"))

fisher_start <- Sys.time()

# Try to use ASE as initialization if possible
# Note: fit_grdpg_cov will run ASE internally anyway, but this documents our intent
fit_fisher <- fit_grdpg_cov(Aselect, B, d = d, p = d, q = 0,
                             maxit = 30, tol = 0.01, tau = tau)

fisher_time <- as.numeric(difftime(Sys.time(), fisher_start, units = "secs"))

cat(sprintf("\n✓ Fisher completed in %.2f seconds (%.2f hours)\n",
            fisher_time, fisher_time/3600))
cat(sprintf("  Converged: %s\n", ifelse(fit_fisher$converged, "YES", "NO")))
cat(sprintf("  Final loss: %.4f\n\n", tail(fit_fisher$loss, 1)))

X_fisher <- fit_fisher$X

# Evaluate Fisher
cat("Evaluating Fisher...\n")
perf_fisher <- evaluate_clustering(X_fisher, labelselect, "FISHER")

# SAVE COMPLETE RESULTS
complete_results <- list(
  data = list(n = n, p = p, k = length(unique(labelselect)), labels = labelselect),
  latent_positions = list(fisher = X_fisher, ase = X_ase, ose = X_ose),
  performance = list(fisher = perf_fisher, ase = perf_ase, ose = perf_ose),
  timing = list(
    fisher = fisher_time,
    ase = ase_time,
    ose = ose_time,
    eval = eval_time
  ),
  fisher_info = list(
    converged = fit_fisher$converged,
    loss = fit_fisher$loss
  )
)

saveRDS(complete_results, "lastfm_complete_results.rds")

cat("\n============================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  ASE:    %.2fs, Sil=%.3f, ARI=%.3f\n",
            ase_time, perf_ase$silhouette, perf_ase$ari))
cat(sprintf("  OSE:    %.2fs, Sil=%.3f, ARI=%.3f\n",
            ose_time, perf_ose$silhouette, perf_ose$ari))
cat(sprintf("  FISHER: %.2fh, Sil=%.3f, ARI=%.3f\n",
            fisher_time/3600, perf_fisher$silhouette, perf_fisher$ari))

cat("\nResults saved:\n")
cat("  - lastfm_ase_ose_results.rds (early save)\n")
cat("  - lastfm_complete_results.rds (final with Fisher)\n")
cat("============================================================================\n")
