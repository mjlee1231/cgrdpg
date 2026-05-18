#!/usr/bin/env Rscript
# Install required R packages for LastFM analysis on HPC

cat("============================================================================\n")
cat("  Installing required R packages for LastFM analysis\n")
cat("============================================================================\n\n")

# List of required packages
packages <- c(
  "igraph",
  "Matrix",
  "irlba",
  "ggplot2",
  "cluster",
  "cgrdpg"  # Your package
)

cat("Packages to install:\n")
cat(paste("  -", packages, collapse="\n"), "\n\n")

# Function to install if not already installed
install_if_needed <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    if (pkg == "cgrdpg") {
      # Install cgrdpg from source
      install.packages("cgrdpg", repos = NULL, type = "source")
    } else {
      # Install from CRAN
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    cat(sprintf("  %s installed successfully!\n\n", pkg))
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# Install each package
for (pkg in packages) {
  tryCatch({
    install_if_needed(pkg)
  }, error = function(e) {
    cat(sprintf("ERROR installing %s: %s\n", pkg, e$message))
  })
}

cat("\n============================================================================\n")
cat("  Installation complete! Testing packages...\n")
cat("============================================================================\n\n")

# Test loading each package
for (pkg in packages) {
  success <- require(pkg, character.only = TRUE, quietly = TRUE)
  status <- if (success) "OK" else "FAILED"
  cat(sprintf("  %-15s %s\n", pkg, status))
}

cat("\n============================================================================\n")
cat("  Done!\n")
cat("============================================================================\n")
