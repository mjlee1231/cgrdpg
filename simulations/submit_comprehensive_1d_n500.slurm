#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=1d_n500_100
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=1d_n500_100reps_%j.out
#SBATCH --error=1d_n500_100reps_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# COMPREHENSIVE 1D SIMULATION: FISHER vs ASE vs OSE (n=500, 100 reps)
# ============================================================================
# Single job running all 100 replications sequentially
# Compares all three methods on the same datasets
# Expected time: ~30-60 minutes for 100 reps
# ============================================================================

# Load R module
module load r/4.3.1

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Simulation: 1D comprehensive (Fisher vs ASE vs OSE)"
echo "Sample size: n=500, 100 replications"
echo "========================================"

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script
Rscript comprehensive_1d_n500_100reps.R

# Check if script succeeded
EXIT_CODE=$?

# Record end time and compute duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"

if [ $EXIT_CODE -eq 0 ]; then
  echo "Status: SUCCESS"
  echo "Results saved to: comprehensive_1d_n500_100reps_results.rds"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
