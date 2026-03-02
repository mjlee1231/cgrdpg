#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=1d_n500
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --array=1-100
#SBATCH --output=1d_n500_array/rep_%03a_%j.out
#SBATCH --error=1d_n500_array/rep_%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# COMPREHENSIVE 1D SIMULATION ARRAY: FISHER vs ASE vs OSE (n=500)
# ============================================================================
# 100 parallel jobs, each running 1 replication
# Each job: ~5 minutes
# Total wall time: ~5 minutes (plus queue wait)
# ============================================================================

# Load R module
module load r/4.3.1

# Create output directory if needed
mkdir -p 1d_n500_array

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Simulation: 1D Fisher vs ASE vs OSE"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "========================================"

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as rep_id
Rscript comprehensive_1d_n500_single.R $SLURM_ARRAY_TASK_ID

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
