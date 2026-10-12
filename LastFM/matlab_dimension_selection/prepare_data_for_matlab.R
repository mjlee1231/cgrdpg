#!/usr/bin/env Rscript
# Convert LastFM data to MATLAB format

library(R.matlab)

cat("Converting LastFM data to MATLAB format...\n")

# Load R data
load("../CovariatesMatrix.Rdata")

# Save as MATLAB .mat file
writeMat("CovariatesMatrix.mat", sparseX = as.matrix(sparseX))

cat("✓ CovariatesMatrix.mat created\n")
cat("✓ Data conversion complete\n")
