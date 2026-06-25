#!/usr/bin/env Rscript
# Dimension Selection for LastFM Data using BIC
# Fits GRDPG with covariates for d = 1, 2, 3, 4, 5
# Computes BIC and selects optimal dimension

library(cgrdpg)
library(ggplot2)

cat("============================================================================\n")
cat("  DIMENSION SELECTION FOR LASTFM DATA USING BIC\n")
cat("  Fitting GRDPG with covariates for d = 1, 2, 3, 4, 5\n")
cat("============================================================================\n\n")

# ============================================================================
# PREPROCESSING (same as analyze_lastfm_modified.R)
# ============================================================================
cat("STEP 1: Loading and preprocessing data...\n")

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

cat(sprintf("  FINAL: n=%d, p=%d, k=%d countries\n\n",
            n, p, length(unique(labelselect))))

# Compute tau
P_temp <- mean(Aselect[upper.tri(Aselect)])
tau <- min(P_temp / 10, 0.001)
cat(sprintf("  Network density: %.4f\n", P_temp))
cat(sprintf("  tau: %.6f\n\n", tau))

# ============================================================================
# DIMENSION SELECTION
# ============================================================================
cat("============================================================================\n")
cat("STEP 2: Fitting models for different dimensions\n")
cat("============================================================================\n\n")

# Dimensions to test
d_values <- 1:5
maxit <- 10
tol <- 0.01

# Storage for results
results <- data.frame(
  d = integer(),
  neg_log_lik = numeric(),
  n_params = integer(),
  BIC = numeric(),
  AIC = numeric(),
  converged = logical(),
  iters = integer(),
  time_sec = numeric()
)

for (d in d_values) {
  cat(sprintf("----------------------------------------------\n"))
  cat(sprintf("Fitting model with d = %d\n", d))
  cat(sprintf("----------------------------------------------\n"))

  start_time <- Sys.time()

  # Fit model
  fit <- fit_grdpg_cov(Aselect, B, d = d, p = d, q = 0,
                       maxit = maxit, tol = tol, tau = tau)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Extract final objective (negative log pseudo-likelihood)
  if (length(fit$history$objective) > 0) {
    neg_log_lik <- tail(fit$history$objective, 1)
  } else {
    neg_log_lik <- NA
  }

  # Compute number of parameters
  # X: n * d parameters (accounting for rotational invariance: subtract d^2)
  # Z: p * d parameters
  # Total: n*d + p*d - d^2
  n_params <- n * d + p * d - d * d

  # Compute BIC: BIC = 2*neg_log_lik + k*log(n)
  # For network data, effective sample size could be n*(n-1)/2 (edges)
  # or n (nodes). We'll use n as it's more conservative.
  BIC <- 2 * neg_log_lik + n_params * log(n)

  # Compute AIC: AIC = 2*neg_log_lik + 2*k
  AIC <- 2 * neg_log_lik + 2 * n_params

  cat(sprintf("  Converged: %s\n", ifelse(fit$converged, "YES", "NO")))
  cat(sprintf("  Iterations: %d\n", fit$iters))
  cat(sprintf("  Time: %.1f sec\n", elapsed))
  cat(sprintf("  Neg-log-lik: %.2f\n", neg_log_lik))
  cat(sprintf("  Num params: %d\n", n_params))
  cat(sprintf("  BIC: %.2f\n", BIC))
  cat(sprintf("  AIC: %.2f\n\n", AIC))

  # Store results
  results <- rbind(results, data.frame(
    d = d,
    neg_log_lik = neg_log_lik,
    n_params = n_params,
    BIC = BIC,
    AIC = AIC,
    converged = fit$converged,
    iters = fit$iters,
    time_sec = elapsed
  ))
}

# ============================================================================
# RESULTS AND VISUALIZATION
# ============================================================================
cat("============================================================================\n")
cat("DIMENSION SELECTION RESULTS\n")
cat("============================================================================\n\n")

print(results, row.names = FALSE)

# Find optimal dimensions
optimal_bic <- results$d[which.min(results$BIC)]
optimal_aic <- results$d[which.min(results$AIC)]

cat(sprintf("\n\nOPTIMAL DIMENSION:\n"))
cat(sprintf("  By BIC: d = %d (BIC = %.2f)\n", optimal_bic, min(results$BIC, na.rm=TRUE)))
cat(sprintf("  By AIC: d = %d (AIC = %.2f)\n", optimal_aic, min(results$AIC, na.rm=TRUE)))

# Save results
saveRDS(results, "dimension_selection_results.rds")
write.csv(results, "dimension_selection_results.csv", row.names = FALSE)
cat("\n\nResults saved to:\n")
cat("  - dimension_selection_results.rds\n")
cat("  - dimension_selection_results.csv\n\n")

# ============================================================================
# PLOTS
# ============================================================================
cat("Creating plots...\n")

# Plot 1: BIC and AIC vs dimension
results_long <- data.frame(
  d = rep(results$d, 2),
  value = c(results$BIC, results$AIC),
  criterion = rep(c("BIC", "AIC"), each = nrow(results))
)

p1 <- ggplot(results_long, aes(x = d, y = value, color = criterion, group = criterion)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = optimal_bic, linetype = "dashed", color = "#F8766D", alpha = 0.7) +
  geom_vline(xintercept = optimal_aic, linetype = "dashed", color = "#00BFC4", alpha = 0.7) +
  scale_x_continuous(breaks = d_values) +
  scale_color_manual(values = c("BIC" = "#F8766D", "AIC" = "#00BFC4")) +
  labs(
    title = "Model Selection: BIC and AIC vs Embedding Dimension",
    subtitle = sprintf("Optimal: BIC selects d=%d, AIC selects d=%d", optimal_bic, optimal_aic),
    x = "Embedding Dimension (d)",
    y = "Information Criterion",
    color = "Criterion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave("dimension_selection_criteria.pdf", p1, width = 8, height = 6)
ggsave("dimension_selection_criteria.png", p1, width = 8, height = 6, dpi = 300)

# Plot 2: Negative log-likelihood vs dimension
p2 <- ggplot(results, aes(x = d, y = neg_log_lik)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = d_values) +
  labs(
    title = "Negative Log Pseudo-Likelihood vs Dimension",
    subtitle = "Lower is better (better fit to data)",
    x = "Embedding Dimension (d)",
    y = "Negative Log Pseudo-Likelihood"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("dimension_selection_loglik.pdf", p2, width = 8, height = 6)
ggsave("dimension_selection_loglik.png", p2, width = 8, height = 6, dpi = 300)

# Plot 3: Number of parameters vs dimension
p3 <- ggplot(results, aes(x = d, y = n_params)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  scale_x_continuous(breaks = d_values) +
  labs(
    title = "Model Complexity vs Dimension",
    subtitle = sprintf("Parameters = n*d + p*d - d² (n=%d, p=%d)", n, p),
    x = "Embedding Dimension (d)",
    y = "Number of Parameters"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("dimension_selection_params.pdf", p3, width = 8, height = 6)
ggsave("dimension_selection_params.png", p3, width = 8, height = 6, dpi = 300)

cat("\nPlots saved:\n")
cat("  - dimension_selection_criteria.pdf/png\n")
cat("  - dimension_selection_loglik.pdf/png\n")
cat("  - dimension_selection_params.pdf/png\n")

cat("\n============================================================================\n")
cat("DIMENSION SELECTION COMPLETE\n")
cat("============================================================================\n")
cat(sprintf("\nRECOMMENDATION: Use d = %d for subsequent analysis\n", optimal_bic))
cat("============================================================================\n")
