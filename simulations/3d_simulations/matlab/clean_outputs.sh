#!/bin/bash
# Clean all output files and directories

echo "Cleaning MATLAB outputs..."

# Remove results directory and recreate empty
rm -rf results/
mkdir -p results

# Remove any stray .mat files
find . -maxdepth 1 -name "*.mat" -delete

# Remove log directories if they exist
rm -rf logs/
rm -rf errors/

echo "Local cleanup complete!"
echo ""
echo "To clean HPC outputs, run on HPC:"
echo "  cd ~/cgrdpg/simulations/3d_simulations/matlab"
echo "  rm -rf results/ logs/ errors/"
echo "  mkdir -p results logs errors"
