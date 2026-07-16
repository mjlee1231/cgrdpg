#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_3d_cov_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=1-100
#SBATCH --output=logs_matlab_3d/rep%03a_%j.out
#SBATCH --error=logs_matlab_3d/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Create output directories
mkdir -p logs_matlab_3d
mkdir -p results_matlab_3d_coverage

# Load MATLAB module
module load matlab/2025a

# Run single replication
matlab -batch "coverage_3d_single_rep(${SLURM_ARRAY_TASK_ID})"
