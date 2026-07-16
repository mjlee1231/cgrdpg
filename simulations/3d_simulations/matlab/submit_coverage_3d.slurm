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

# Run all 100 replications in a loop
matlab -batch "
fprintf('Running 3D Coverage: 100 Replications\n');
fprintf('======================================\n\n');
addpath('core');
total_start = tic;
for rep = 1:100
    fprintf('\n--- Rep %d/100 ---\n', rep);
    try
        coverage_3d_single_rep(rep);
    catch ME
        fprintf('ERROR in rep %d: %s\n', rep, ME.message);
    end
end
fprintf('\nTotal time: %.1f min\n', toc(total_start)/60);
fprintf('Now run: matlab -batch \"aggregate_coverage_3d\"\n');
"
