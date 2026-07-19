#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=test_latent_design
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs_test_design/test_%j.out
#SBATCH --error=logs_test_design/test_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Create output directory
mkdir -p logs_test_design

# Load MATLAB module
module load matlab/2025a

# Run test comparing old vs new latent position designs
matlab -batch "test_new_latent_design"
