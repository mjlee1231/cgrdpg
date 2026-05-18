#!/usr/bin/env Rscript
# Visualize LastFM GRDPG analysis results
# Creates plots comparing Fisher vs ASE vs OSE latent positions

library(ggplot2)
library(gridExtra)
library(reshape2)

cat("============================================================================\n")
cat("  VISUALIZING LASTFM GRDPG RESULTS\n")
cat("============================================================================\n\n")

# Load results
results <- readRDS("lastfm_grdpg_analysis_results.rds")

X_fisher <- results$latent_positions$fisher
X_ase <- results$latent_positions$ase
X_ose <- results$latent_positions$ose
labels <- results$data$labels

n <- nrow(X_fisher)
d <- ncol(X_fisher)

cat(sprintf("Loaded results: n=%d users, d=%d dimensions\n", n, d))
cat(sprintf("Countries: %s\n\n", paste(unique(labels), collapse=", ")))

# Color palette for countries
prettyColors <- c("turquoise4", "red", "cyan", "deeppink4", "violet", "blue", "orange", "green")
label_unique <- unique(labels)
k <- length(label_unique)

# ============================================================================
# 2D SCATTER PLOTS (First 2 dimensions)
# ============================================================================
cat("Creating 2D scatter plots...\n")

# Prepare data frames
df_fisher <- data.frame(
  Dim1 = X_fisher[, 1],
  Dim2 = X_fisher[, 2],
  Country = as.factor(labels),
  Method = "Fisher"
)

df_ase <- data.frame(
  Dim1 = X_ase[, 1],
  Dim2 = X_ase[, 2],
  Country = as.factor(labels),
  Method = "ASE"
)

df_ose <- data.frame(
  Dim1 = X_ose[, 1],
  Dim2 = X_ose[, 2],
  Country = as.factor(labels),
  Method = "OSE"
)

# Individual plots
p1 <- ggplot(df_fisher, aes(x = Dim1, y = Dim2, color = Country)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = prettyColors[1:k]) +
  labs(title = "FISHER (with covariates)",
       x = "Dimension 1", y = "Dimension 2") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(df_ase, aes(x = Dim1, y = Dim2, color = Country)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = prettyColors[1:k]) +
  labs(title = "ASE (no covariates)",
       x = "Dimension 1", y = "Dimension 2") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

p3 <- ggplot(df_ose, aes(x = Dim1, y = Dim2, color = Country)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = prettyColors[1:k]) +
  labs(title = "OSE (no covariates)",
       x = "Dimension 1", y = "Dimension 2") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save individual plots
ggsave("lastfm_fisher_2d.pdf", p1, width = 8, height = 6)
ggsave("lastfm_ase_2d.pdf", p2, width = 8, height = 6)
ggsave("lastfm_ose_2d.pdf", p3, width = 8, height = 6)

# Combined plot
pdf("lastfm_comparison_2d.pdf", width = 18, height = 6)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

cat("  Saved: lastfm_fisher_2d.pdf, lastfm_ase_2d.pdf, lastfm_ose_2d.pdf\n")
cat("  Saved: lastfm_comparison_2d.pdf\n\n")

# ============================================================================
# 3D SCATTER PLOTS (if d >= 3)
# ============================================================================
if (d >= 3) {
  cat("Creating 3D scatter plots...\n")

  # Pairs plot for first 3 dimensions
  df_fisher_3d <- data.frame(
    Dim1 = X_fisher[, 1],
    Dim2 = X_fisher[, 2],
    Dim3 = X_fisher[, 3],
    Country = as.factor(labels)
  )

  df_ase_3d <- data.frame(
    Dim1 = X_ase[, 1],
    Dim2 = X_ase[, 2],
    Dim3 = X_ase[, 3],
    Country = as.factor(labels)
  )

  df_ose_3d <- data.frame(
    Dim1 = X_ose[, 1],
    Dim2 = X_ose[, 2],
    Dim3 = X_ose[, 3],
    Country = as.factor(labels)
  )

  # Pairs plots
  pdf("lastfm_fisher_pairs.pdf", width = 10, height = 10)
  pairs(df_fisher_3d[, 1:3],
        col = prettyColors[as.numeric(as.factor(labels))],
        pch = 16, main = "FISHER (with covariates) - First 3 Dimensions")
  dev.off()

  pdf("lastfm_ase_pairs.pdf", width = 10, height = 10)
  pairs(df_ase_3d[, 1:3],
        col = prettyColors[as.numeric(as.factor(labels))],
        pch = 16, main = "ASE (no covariates) - First 3 Dimensions")
  dev.off()

  pdf("lastfm_ose_pairs.pdf", width = 10, height = 10)
  pairs(df_ose_3d[, 1:3],
        col = prettyColors[as.numeric(as.factor(labels))],
        pch = 16, main = "OSE (no covariates) - First 3 Dimensions")
  dev.off()

  cat("  Saved: lastfm_fisher_pairs.pdf, lastfm_ase_pairs.pdf, lastfm_ose_pairs.pdf\n\n")
}

# ============================================================================
# PERFORMANCE COMPARISON BARPLOTS
# ============================================================================
cat("Creating performance comparison plots...\n")

perf_df <- data.frame(
  Method = rep(c("FISHER", "ASE", "OSE"), 3),
  Metric = rep(c("Silhouette", "ARI", "BSS/TSS"), each = 3),
  Value = c(
    results$performance$fisher$silhouette,
    results$performance$ase$silhouette,
    results$performance$ose$silhouette,
    results$performance$fisher$ari,
    results$performance$ase$ari,
    results$performance$ose$ari,
    results$performance$fisher$bss_tss,
    results$performance$ase$bss_tss,
    results$performance$ose$bss_tss
  )
)

p_perf <- ggplot(perf_df, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c("FISHER" = "darkgreen", "ASE" = "steelblue", "OSE" = "darkorange")) +
  labs(title = "Performance Comparison: Country Label Alignment",
       y = "Score (higher is better)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("lastfm_performance_comparison.pdf", p_perf, width = 10, height = 4)
cat("  Saved: lastfm_performance_comparison.pdf\n\n")

# ============================================================================
# TIMING COMPARISON
# ============================================================================
cat("Creating timing comparison plot...\n")

timing_df <- data.frame(
  Method = c("FISHER", "ASE", "OSE"),
  Time = c(
    results$timing$fisher,
    results$timing$ase,
    results$timing$ose
  )
)

p_timing <- ggplot(timing_df, aes(x = Method, y = Time, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FISHER" = "darkgreen", "ASE" = "steelblue", "OSE" = "darkorange")) +
  labs(title = "Computational Time Comparison",
       y = "Time (seconds)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("lastfm_timing_comparison.pdf", p_timing, width = 6, height = 4)
cat("  Saved: lastfm_timing_comparison.pdf\n\n")

# ============================================================================
# SUMMARY REPORT
# ============================================================================
cat("============================================================================\n")
cat("  VISUALIZATION COMPLETE\n")
cat("============================================================================\n\n")

cat("Generated files:\n")
cat("  - lastfm_fisher_2d.pdf\n")
cat("  - lastfm_ase_2d.pdf\n")
cat("  - lastfm_ose_2d.pdf\n")
cat("  - lastfm_comparison_2d.pdf\n")
if (d >= 3) {
  cat("  - lastfm_fisher_pairs.pdf\n")
  cat("  - lastfm_ase_pairs.pdf\n")
  cat("  - lastfm_ose_pairs.pdf\n")
}
cat("  - lastfm_performance_comparison.pdf\n")
cat("  - lastfm_timing_comparison.pdf\n\n")

cat("QUICK SUMMARY:\n")
cat(sprintf("  Best Silhouette: %s (%.4f)\n",
            c("FISHER", "ASE", "OSE")[which.max(c(results$performance$fisher$silhouette,
                                                    results$performance$ase$silhouette,
                                                    results$performance$ose$silhouette))],
            max(results$performance$fisher$silhouette,
                results$performance$ase$silhouette,
                results$performance$ose$silhouette)))

cat(sprintf("  Best ARI: %s (%.4f)\n",
            c("FISHER", "ASE", "OSE")[which.max(c(results$performance$fisher$ari,
                                                    results$performance$ase$ari,
                                                    results$performance$ose$ari))],
            max(results$performance$fisher$ari,
                results$performance$ase$ari,
                results$performance$ose$ari)))

cat(sprintf("  Fastest: %s (%.2f sec)\n",
            c("FISHER", "ASE", "OSE")[which.min(c(results$timing$fisher,
                                                    results$timing$ase,
                                                    results$timing$ose))],
            min(results$timing$fisher,
                results$timing$ase,
                results$timing$ose)))

cat("\n============================================================================\n")
