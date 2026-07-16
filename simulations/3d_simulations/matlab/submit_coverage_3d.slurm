#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_3d_cov
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=matlab_3d_coverage_%j.out
#SBATCH --error=matlab_3d_coverage_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Create output directory
mkdir -p results_matlab_3d_coverage

# Load MATLAB module
module load matlab/2025a

# Run all 100 replications
matlab -batch "run_all_3d"
