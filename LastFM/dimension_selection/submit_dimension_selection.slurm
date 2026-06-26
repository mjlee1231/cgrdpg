#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=dimsel_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --array=1-5
#SBATCH --output=logs_tau0.001/dimension_d%a_%j.out
#SBATCH --error=errors_tau0.001/dimension_d%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# Dimension Selection for LastFM Data - PARALLEL EXECUTION (tau=0.001)
# ============================================================================
# Runs fit_grdpg_cov_parallel for d = 1, 2, 3, 4, 5 in parallel
# Each array job handles one dimension with 32 cores
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
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Dimension: d = $SLURM_ARRAY_TASK_ID"
echo "Resources: 32 CPUs, 128GB RAM, 4 hours"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

# Run R script with dimension as argument
Rscript dimension_selection_single.R $SLURM_ARRAY_TASK_ID

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
