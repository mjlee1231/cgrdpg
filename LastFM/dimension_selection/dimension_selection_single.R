#!/usr/bin/env Rscript
# Dimension Selection for LastFM - SINGLE DIMENSION (for HPC parallel execution)
# Takes dimension as command-line argument

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript dimension_selection_single.R <dimension>")
}

d <- as.integer(args[1])
if (is.na(d) || d < 1 || d > 5) {
  stop("Dimension must be between 1 and 5")
}

library(cgrdpg)

cat("============================================================================\n")
cat(sprintf("  DIMENSION SELECTION: d = %d\n", d))
cat("  HPC Parallel Execution\n")
cat("============================================================================\n\n")

cat(sprintf("Job started: %s\n", Sys.time()))
cat(sprintf("Dimension: d = %d\n\n", d))

# ============================================================================
# PREPROCESSING
# ============================================================================
cat("STEP 1: Loading and preprocessing data...\n")
flush.console()

load("../CovariatesMatrix.Rdata")
X <- as.matrix(sparseX)
edgelist <- read.csv("../lastfm_asia_edges.csv") + 1
A <- matrix(0, nrow(X), nrow(X))
for (i in 1:nrow(edgelist)) A[edgelist[i,1], edgelist[i,2]] <- 1
A <- A + t(A)

label <- read.csv('../lastfm_asia_target.csv')$target + 1
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

n <- nrow(Aselect)
p <- ncol(Xselect)
B <- t(Xselect)

cat(sprintf("  n=%d, p=%d, k=%d countries\n", n, p, length(unique(labelselect))))
flush.console()

# Compute tau
P_temp <- mean(Aselect[upper.tri(Aselect)])
tau <- min(P_temp / 10, 0.001)
cat(sprintf("  Network density: %.4f, tau: %.6f\n\n", P_temp, tau))
flush.console()

# ============================================================================
# FIT MODEL FOR THIS DIMENSION
# ============================================================================
cat("============================================================================\n")
cat(sprintf("STEP 2: Fitting model with d = %d\n", d))
cat("============================================================================\n\n")
flush.console()

maxit <- 10
tol <- 0.01

# Get number of cores from SLURM
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (ncores <= 1) ncores <- max(1, parallel::detectCores() - 1)

cat(sprintf("Starting Fisher-scoring with parallel computation...\n"))
cat(sprintf("  maxit = %d, tol = %.3f, ncores = %d\n", maxit, tol, ncores))
cat(sprintf("  Fit started: %s\n\n", Sys.time()))
flush.console()

start_time <- Sys.time()

# Fit model with parallel computation
fit <- tryCatch(
  fit_grdpg_cov_parallel(Aselect, B, d = d, p = d, q = 0,
                         maxit = maxit, tol = tol, tau = tau, ncores = ncores),
  error = function(e) {
    cat("Parallel failed, falling back to sequential...\n")
    flush.console()
    fit_grdpg_cov(Aselect, B, d = d, p = d, q = 0,
                  maxit = maxit, tol = tol, tau = tau)
  }
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n✓ Fitting completed in %.2f seconds (%.2f hours)\n", elapsed, elapsed/3600))
cat(sprintf("  Converged: %s\n", ifelse(fit$converged, "YES", "NO")))
cat(sprintf("  Iterations: %d\n", fit$iters))
flush.console()

# Extract final objective (negative log pseudo-likelihood)
if (length(fit$history$objective) > 0) {
  neg_log_lik <- tail(fit$history$objective, 1)
  cat(sprintf("  Final neg-log-lik: %.4f\n\n", neg_log_lik))

  cat("  Loss trajectory:\n")
  for (i in 1:length(fit$history$objective)) {
    cat(sprintf("    Iter %d: %.6f\n", i-1, fit$history$objective[i]))
  }
  cat("\n")
  flush.console()
} else {
  neg_log_lik <- NA
  cat("  WARNING: No objective history recorded\n\n")
}

if (length(fit$history$max_row_change) > 0) {
  cat("  Max row change trajectory:\n")
  for (i in 1:length(fit$history$max_row_change)) {
    cat(sprintf("    Iter %d: %.6f\n", i, fit$history$max_row_change[i]))
  }
  cat("\n")
  flush.console()
}

# ============================================================================
# COMPUTE INFORMATION CRITERIA
# ============================================================================
cat("============================================================================\n")
cat("STEP 3: Computing information criteria\n")
cat("============================================================================\n\n")

# Compute number of parameters
# X: n * d parameters (accounting for rotational invariance: subtract d^2)
# Z: p * d parameters
# Total: n*d + p*d - d^2
n_params <- n * d + p * d - d * d

# Compute BIC: BIC = 2*neg_log_lik + k*log(n)
BIC <- 2 * neg_log_lik + n_params * log(n)

# Compute AIC: AIC = 2*neg_log_lik + 2*k
AIC <- 2 * neg_log_lik + 2 * n_params

cat(sprintf("  Neg-log-likelihood: %.4f\n", neg_log_lik))
cat(sprintf("  Number of parameters: %d\n", n_params))
cat(sprintf("  BIC: %.4f\n", BIC))
cat(sprintf("  AIC: %.4f\n\n", AIC))
flush.console()

# ============================================================================
# SAVE RESULTS
# ============================================================================
cat("============================================================================\n")
cat("STEP 4: Saving results\n")
cat("============================================================================\n\n")

result <- list(
  d = d,
  neg_log_lik = neg_log_lik,
  n_params = n_params,
  BIC = BIC,
  AIC = AIC,
  converged = fit$converged,
  iters = fit$iters,
  time_sec = elapsed,
  data_info = list(n = n, p = p, k = length(unique(labelselect))),
  fit_history = list(
    objective = fit$history$objective,
    max_row_change = fit$history$max_row_change
  ),
  latent_positions = fit$X,
  covariate_effects = fit$Z,
  timestamp = Sys.time()
)

output_file <- sprintf("results_d%d.rds", d)
saveRDS(result, output_file)

cat(sprintf("✓ Results saved to: %s\n\n", output_file))
flush.console()

# ============================================================================
# SUMMARY
# ============================================================================
cat("============================================================================\n")
cat("JOB COMPLETE\n")
cat("============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  Dimension: d = %d\n", d))
cat(sprintf("  Converged: %s (%d iterations)\n", ifelse(fit$converged, "YES", "NO"), fit$iters))
cat(sprintf("  Time: %.2f hours\n", elapsed/3600))
cat(sprintf("  Neg-log-lik: %.4f\n", neg_log_lik))
cat(sprintf("  BIC: %.4f\n", BIC))
cat(sprintf("  AIC: %.4f\n", AIC))
cat(sprintf("  Output: %s\n", output_file))
cat(sprintf("\nJob completed: %s\n", Sys.time()))
cat("============================================================================\n")
flush.console()
