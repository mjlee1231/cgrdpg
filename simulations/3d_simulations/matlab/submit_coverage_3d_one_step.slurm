#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_one_step_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=1-100
#SBATCH --output=logs_matlab_one_step/rep%03a_%j.out
#SBATCH --error=logs_matlab_one_step/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# One-step cgrdpg estimators (Batch vs Jacobi)
# Should be MUCH faster than full optimization (~30 sec vs 2 min)

# Create output directories
mkdir -p logs_matlab_one_step
mkdir -p results_matlab_3d_one_step

# Load MATLAB module
module load matlab/2025a

# Run single replication
matlab -batch "coverage_3d_single_rep_one_step(${SLURM_ARRAY_TASK_ID})"
