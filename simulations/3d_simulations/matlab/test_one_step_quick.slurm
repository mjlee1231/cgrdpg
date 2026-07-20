#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=test_one_step
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=test_one_step_output.txt
#SBATCH --error=test_one_step_error.txt

# Quick test to verify one-step bug fix

# Load MATLAB module
module load matlab/2025a

# Run debug script
matlab -batch "test_one_step_debug"
