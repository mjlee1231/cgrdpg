#!/usr/bin/env Rscript
# Install required R packages for LastFM analysis on HPC

cat("============================================================================\n")
cat("Installing R packages for LastFM analysis\n")
cat("============================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of required packages
required_packages <- c(
  "igraph",
  "Matrix",
  "irlba",
  "ggplot2",
  "cluster",
  "foreach",
  "doParallel"
)

cat("Required packages:\n")
for (pkg in required_packages) {
  cat(sprintf("  - %s\n", pkg))
}
cat("\n")

# Install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("  ✓ %s installed successfully\n", pkg))
    } else {
      cat(sprintf("  ✗ Failed to install %s\n", pkg))
    }
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n============================================================================\n")
cat("Checking cgrdpg package...\n")
cat("============================================================================\n\n")

# Check if cgrdpg is installed
if (require("cgrdpg", quietly = TRUE)) {
  cat("✓ cgrdpg package is installed\n")
  cat(sprintf("  Version: %s\n", packageVersion("cgrdpg")))
} else {
  cat("✗ cgrdpg package is NOT installed\n")
  cat("\nPlease install cgrdpg manually:\n")
  cat("  1. cd ~/cgrdpg\n")
  cat("  2. module load r/4.4.1\n")
  cat("  3. R CMD INSTALL .\n")
}

cat("\n============================================================================\n")
cat("Installation complete!\n")
cat("============================================================================\n")
