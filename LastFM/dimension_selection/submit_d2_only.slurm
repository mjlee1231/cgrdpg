#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=dimsel_d2
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs_tau0.001/dimension_d2_rerun_%j.out
#SBATCH --error=errors_tau0.001/dimension_d2_rerun_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# Dimension Selection - D=2 ONLY (Extended time for slow dimension)
# ============================================================================
# d=2 is the slowest dimension, needs 8 hours instead of 4
# ============================================================================

module load r/4.3.1

# Create directories for outputs (tau=0.001)
mkdir -p logs_tau0.001
mkdir -p errors_tau0.001
mkdir -p results_tau0.001

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Dimension: d = 2 (rerun with 8-hour limit)"
echo "Resources: 32 CPUs, 128GB RAM, 8 hours"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

# Run R script with d=2
Rscript dimension_selection_single.R 2

EXIT_CODE=$?

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "Runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"

if [ $EXIT_CODE -eq 0 ]; then
  echo "Status: SUCCESS"
  echo "Dimension d=2 complete"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
