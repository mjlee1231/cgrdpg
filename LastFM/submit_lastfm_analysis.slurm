#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=lastfm_grdpg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=lastfm_grdpg_%j.out
#SBATCH --error=lastfm_grdpg_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# LASTFM ASIA SOCIAL NETWORK ANALYSIS WITH GRDPG
# ============================================================================
# Comparing Fisher-scoring (with covariates) vs ASE vs OSE (no covariates)
# Expected runtime: ~30-60 minutes
# ============================================================================

# Load R module
module load r/4.3.1

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32GB"
echo "Analysis: LastFM GRDPG (Fisher vs ASE vs OSE)"
echo "========================================"

# Change to LastFM directory
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script
Rscript analyze_lastfm_grdpg_hpc.R

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
  echo "Results saved to: lastfm_grdpg_analysis_results_hpc.rds"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
