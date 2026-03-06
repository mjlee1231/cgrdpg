#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=lastfm_modified
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=lastfm_modified_%j.out
#SBATCH --error=lastfm_modified_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# LASTFM ANALYSIS WITH MODIFIED FISHER-SCORING (maxit=10)
# ============================================================================
# Modification: Clipped Fisher information for better conditioning
# Expected runtime: ~2-4 hours
# ============================================================================

# Load R module
module load r/4.3.1

# Set environment variables to enable multi-threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLIS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 200GB"
echo "Time limit: 12 hours"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "Analysis: LastFM with MODIFIED Fisher-scoring (maxit=10)"
echo "Modification: Clipped Fisher info for better conditioning"
echo "========================================"

# Change to LastFM directory
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script
Rscript analyze_lastfm_modified_hpc.R

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
  echo "Results saved to:"
  echo "  - lastfm_modified_ase_ose_results_hpc.rds"
  echo "  - lastfm_modified_complete_results_hpc.rds"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
