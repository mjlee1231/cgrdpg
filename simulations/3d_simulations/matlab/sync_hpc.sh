#!/bin/bash
# Sync HPC MATLAB folder with clean laptop version
# Run this on HPC: bash sync_hpc.sh

echo "========================================="
echo "Syncing HPC MATLAB folder from GitHub"
echo "========================================="

# 1. Pull latest changes from GitHub
echo -e "\n1. Pulling latest changes..."
git fetch origin
git checkout claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN
git pull origin claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN

# 2. Clean up old output/log files
echo -e "\n2. Cleaning up old log and output files..."
rm -f matlab_3d_*.err matlab_3d_*.out
rm -f matlab_2d_*.err matlab_2d_*.out
rm -f matlab_3d_ase_ose_*.err matlab_3d_ase_ose_*.out

# 3. Remove old log directories (keeping the latest results)
echo -e "\n3. Cleaning up old log directories..."
rm -rf errors errors_matlab_2d errors_matlab_3d
rm -rf logs logs_matlab_2d logs_matlab_3d

# 4. List what's left
echo -e "\n4. Current directory contents:"
ls -lh

echo -e "\n========================================="
echo "Sync complete!"
echo "========================================="
echo ""
echo "Files kept:"
echo "  - Core functions (core/*.m)"
echo "  - Main simulation: coverage_3d_single_rep.m"
echo "  - Aggregation: aggregate_coverage_3d.m"
echo "  - SLURM script: submit_coverage_3d.slurm"
echo "  - Results directories (if you want to keep old results)"
echo ""
echo "Files removed:"
echo "  - Old run_all_*.m files"
echo "  - Old SLURM output/error files"
echo "  - coverage_3d_ase_ose_single_rep.m (merged into coverage_3d_single_rep.m)"
echo "  - submit_coverage_3d_ase_ose.slurm (redundant)"
echo ""
echo "To clean old results and start fresh, run:"
echo "  rm -rf results_matlab_3d_ase_ose_cgrdpg"
echo "  rm -rf results_matlab_3d_coverage  # Only if you want to re-run all 100 reps"
echo ""
echo "Ready to submit new jobs:"
echo "  sbatch submit_coverage_3d.slurm"
echo "========================================="
