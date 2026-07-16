#!/bin/bash
#SBATCH --job-name=matlab_2d_cov
#SBATCH --output=logs_matlab_2d/coverage_%a.out
#SBATCH --error=errors_matlab_2d/coverage_%a.err
#SBATCH --array=1-100
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Create output directories
mkdir -p logs_matlab_2d
mkdir -p errors_matlab_2d
mkdir -p results_matlab_2d_coverage

# Load MATLAB module
module load matlab/R2024b

# Run single replication
matlab -batch "coverage_2d_single_rep(${SLURM_ARRAY_TASK_ID})"
