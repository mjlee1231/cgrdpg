#!/usr/bin/env Rscript
# Debug script to inspect .mat file structure

library(R.matlab)

file_path <- "results_matlab_maxit50_tol0.1/results_d1.mat"

if (file.exists(file_path)) {
  cat("Reading:", file_path, "\n\n")
  mat_data <- readMat(file_path)

  cat("Structure of mat_data:\n")
  str(mat_data)

  cat("\n\nNames in mat_data:\n")
  print(names(mat_data))

  cat("\n\nStructure of result:\n")
  str(mat_data$result)

  cat("\n\nDimnames of result:\n")
  print(dimnames(mat_data$result))

  cat("\n\nTrying to extract values:\n")
  result <- mat_data$result

  # Try different indexing approaches
  cat("Number of elements:", length(result), "\n")

  for (i in 1:min(10, length(result))) {
    cat(sprintf("result[[%d]]: ", i))
    print(result[[i]])
  }
} else {
  cat("File not found:", file_path, "\n")
}
