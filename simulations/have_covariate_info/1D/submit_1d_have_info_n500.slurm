#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=1d_have_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=1-100
#SBATCH --output=logs_1d_have_info_n500/rep%03a_%j.out
#SBATCH --error=errors_1d_have_info_n500/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# 1D SIMULATION WITH COVARIATE INFO: FISHER vs ASE vs OSE (n=500)
# ============================================================================
# Z0 = rnorm (covariates contain signal)
# 100 replications, maxit=15, tol=0.05, 16 cores parallel
# Comparing WITH signal vs no_covariate_info (NO signal)
# ============================================================================

# Load R module
module load r/4.3.1

# Create separate directories for results, logs, and errors
mkdir -p results_1d_have_info_n500
mkdir -p logs_1d_have_info_n500
mkdir -p errors_1d_have_info_n500

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scenario: Z0 = rnorm (WITH covariate signal)"
echo "Methods: FISHER-TRUE/PLUGIN, ASE-TRUE/PLUGIN, OSE-TRUE/PLUGIN"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: 16 CPUs, 64GB RAM, maxit=15, tol=0.05"
echo "========================================"

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as rep_id
Rscript vertex_wise_coverage_1d_have_info_n500.R $SLURM_ARRAY_TASK_ID

# Check if script succeeded
EXIT_CODE=$?

# Record end time and compute duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "Runtime: ${MINUTES}m ${SECONDS}s"

if [ $EXIT_CODE -eq 0 ]; then
  echo "Status: SUCCESS"
  echo "Replication $SLURM_ARRAY_TASK_ID complete"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
