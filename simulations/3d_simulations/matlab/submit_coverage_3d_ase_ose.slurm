#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_3d_ase_ose
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=matlab_3d_ase_ose_%j.out
#SBATCH --error=matlab_3d_ase_ose_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Create output directory
mkdir -p results_matlab_3d_ase_ose_cgrdpg

# Load MATLAB module
module load matlab/2025a

# Run all 100 replications
matlab -batch "run_all_3d_ase_ose"
