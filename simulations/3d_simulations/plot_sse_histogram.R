#!/usr/bin/env Rscript
# Plot histogram of SSEs from simulation results

library(ggplot2)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_sse_histogram.R <results_folder_or_file> [output_prefix]\n")
  cat("Example: Rscript plot_sse_histogram.R results_3d_ase_ose_cgrdpg_n1000 sse_n1000\n")
  cat("Example: Rscript plot_sse_histogram.R results_file.rds sse_n1000\n")
  quit(status = 1)
}

input_path <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "sse_histogram"

# Check if input is a directory or file
if (dir.exists(input_path)) {
  # It's a directory - read all .rds files
  cat(sprintf("Reading results from directory: %s\n", input_path))
  rds_files <- list.files(input_path, pattern = "\\.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    cat(sprintf("Error: No .rds files found in directory '%s'\n", input_path))
    quit(status = 1)
  }

  cat(sprintf("Found %d .rds files\n", length(rds_files)))

  # Read all files and combine
  all_results <- lapply(rds_files, function(f) {
    cat(sprintf("  Reading: %s\n", basename(f)))
    readRDS(f)
  })

  # Combine into single structure
  results <- all_results

} else if (file.exists(input_path)) {
  # It's a file
  cat(sprintf("Loading results from file: %s\n", input_path))
  results <- readRDS(input_path)
} else {
  cat(sprintf("Error: Path '%s' not found\n", input_path))
  quit(status = 1)
}

cat(sprintf("Results structure: %s\n", paste(names(results), collapse=", ")))

# Extract SSE data
# Assuming structure: results$sse_ase, results$sse_ose, results$sse_cgrdpg
# Or results might be a list of replications

# Flatten if we read multiple files
if (is.list(results) && length(results) > 0 && is.list(results[[1]]) &&
    !"sse_ase" %in% names(results) && "sse_ase" %in% names(results[[1]])) {
  # Results is a list of lists - flatten it
  cat("Flattening results from multiple files...\n")
  results <- unlist(results, recursive = FALSE)
}

# Try to infer structure
if ("sse_ase" %in% names(results)) {
  # Structure: direct SSE fields (separate)
  cat("Structure: Single result with separate SSE fields\n")
  sse_data <- data.frame(
    SSE = c(results$sse_ase, results$sse_ose, results$sse_cgrdpg),
    Method = rep(c("ASE", "OSE", "CGRDPG"),
                 each = length(results$sse_ase))
  )
} else if ("sse" %in% names(results) && is.numeric(results$sse)) {
  # Structure: single result with named SSE vector
  cat("Structure: Single result with named SSE vector\n")
  sse_vec <- results$sse
  sse_data <- data.frame(
    SSE = as.numeric(sse_vec),
    Method = toupper(names(sse_vec))
  )
} else if (is.list(results) && length(results) > 0) {
  # Structure: list of replications
  cat(sprintf("Structure: List of %d replications\n", length(results)))

  # Check first element structure
  first_elem <- results[[1]]

  if ("sse_ase" %in% names(first_elem)) {
    # Separate SSE fields
    n_reps <- length(results)
    sse_ase <- sapply(results, function(x) {
      if ("sse_ase" %in% names(x)) x$sse_ase else NA
    })
    sse_ose <- sapply(results, function(x) {
      if ("sse_ose" %in% names(x)) x$sse_ose else NA
    })
    sse_cgrdpg <- sapply(results, function(x) {
      if ("sse_cgrdpg" %in% names(x)) x$sse_cgrdpg else NA
    })

    # Remove NAs
    sse_ase <- sse_ase[!is.na(sse_ase)]
    sse_ose <- sse_ose[!is.na(sse_ose)]
    sse_cgrdpg <- sse_cgrdpg[!is.na(sse_cgrdpg)]

    sse_data <- data.frame(
      SSE = c(sse_ase, sse_ose, sse_cgrdpg),
      Method = rep(c("ASE", "OSE", "CGRDPG"),
                   times = c(length(sse_ase), length(sse_ose), length(sse_cgrdpg)))
    )
  } else if ("sse" %in% names(first_elem) && is.numeric(first_elem$sse)) {
    # Named SSE vector structure
    cat("Extracting from named SSE vectors...\n")

    # Extract SSE vectors from each replication
    all_sse <- lapply(results, function(x) {
      if ("sse" %in% names(x)) {
        data.frame(
          SSE = as.numeric(x$sse),
          Method = toupper(names(x$sse))
        )
      } else {
        NULL
      }
    })

    # Remove NULLs and combine
    all_sse <- all_sse[!sapply(all_sse, is.null)]
    sse_data <- do.call(rbind, all_sse)
  } else {
    cat("Error: Cannot identify SSE data structure in list elements\n")
    cat("First element structure:\n")
    print(str(first_elem, max.level = 2))
    quit(status = 1)
  }
} else {
  cat("Error: Cannot identify SSE data structure\n")
  cat("Available fields:\n")
  print(str(results, max.level = 2))
  quit(status = 1)
}

cat(sprintf("\nLoaded %d SSE values across %d methods\n",
            nrow(sse_data), length(unique(sse_data$Method))))

# Summary statistics
cat("\nSummary statistics by method:\n")
print(aggregate(SSE ~ Method, data = sse_data,
                FUN = function(x) c(
                  n = length(x),
                  mean = mean(x),
                  sd = sd(x),
                  median = median(x),
                  min = min(x),
                  max = max(x)
                )))

# Create histogram
cat("\nCreating histogram...\n")

# Use Cairo for HPC compatibility
options(bitmapType = 'cairo')

p <- ggplot(sse_data, aes(x = SSE, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of Sum of Squared Errors (SSE)",
    subtitle = sprintf("Based on %s", basename(input_path)),
    x = "SSE",
    y = "Frequency"
  ) +
  theme(legend.position = "none")

# Save plot
output_file <- sprintf("%s.png", output_prefix)
ggsave(output_file, p, width = 8, height = 10, dpi = 300)
cat(sprintf("✓ Histogram saved to: %s\n", output_file))

# Also create overlaid histogram
p2 <- ggplot(sse_data, aes(x = SSE, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() +
  labs(
    title = "Distribution of Sum of Squared Errors (SSE) - Overlaid",
    subtitle = sprintf("Based on %s", basename(input_path)),
    x = "SSE",
    y = "Frequency"
  ) +
  scale_fill_brewer(palette = "Set2")

output_file2 <- sprintf("%s_overlaid.png", output_prefix)
ggsave(output_file2, p2, width = 10, height = 6, dpi = 300)
cat(sprintf("✓ Overlaid histogram saved to: %s\n", output_file2))

cat("\nDone!\n")
