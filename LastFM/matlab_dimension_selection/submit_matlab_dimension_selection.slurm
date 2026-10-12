#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_dimsel_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-5
#SBATCH --output=logs_matlab/dimension_d%a_%j.out
#SBATCH --error=errors_matlab/dimension_d%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# Dimension Selection for LastFM Data - MATLAB VERSION
# ============================================================================
# Runs MATLAB fit_grdpg_fminunc_surrogate for d = 1, 2, 3, 4, 5 in parallel
# Each array job handles one dimension
# Expected to be faster than R version due to MATLAB's optimized fminunc
# ============================================================================

module load matlab

# Create directories for outputs
mkdir -p logs_matlab
mkdir -p errors_matlab
mkdir -p results_matlab_maxit30_tol0.1

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Dimension: d = $SLURM_ARRAY_TASK_ID"
echo "Resources: 1 CPU, 64GB RAM, 8 hours"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

# Run MATLAB script with dimension as argument
matlab -nodisplay -nosplash -r "dimension_selection_single($SLURM_ARRAY_TASK_ID); exit"

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
  echo "Dimension d=$SLURM_ARRAY_TASK_ID complete"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
