#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_3d_n2000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-100
#SBATCH --output=logs_matlab_3d_n2000/rep%03a_%j.out
#SBATCH --error=logs_matlab_3d_n2000/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Larger sample size: n=2000, p_cov=1000
# Increased memory (32G) and time (4h) for larger problem

# Create output directories
mkdir -p logs_matlab_3d_n2000
mkdir -p results_matlab_3d_coverage_n2000

# Load MATLAB module
module load matlab/2025a

# Run single replication
matlab -batch "coverage_3d_single_rep_n2000(${SLURM_ARRAY_TASK_ID})"
