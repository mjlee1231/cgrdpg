#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=lastfm_d4
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=08:00:00
#SBATCH --output=logs/lastfm_d4_%j.out
#SBATCH --error=errors/lastfm_d4_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# LastFM Analysis: d=4 (2nd optimal dimension by BIC)
# ============================================================================
# Parameters: maxit=20, tol=0.05, tau=0.001
# Using 48 cores for parallel Fisher scoring
# ============================================================================

module load r/4.3.1

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Dimension: d = 4"
echo "Resources: 48 CPUs, 192GB RAM, 8 hours"
echo "Parameters: maxit=20, tol=0.05, tau=0.001"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

# Run R script with d=4
Rscript ../analyze_lastfm_single_dimension.R 4

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
  echo "Analysis for d=4 complete"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
