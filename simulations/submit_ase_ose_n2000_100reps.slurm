#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=ase_n2000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --array=1-100
#SBATCH --output=ase_ose_logs_n2000/ase_ose_n2000_rep%03a_%j.out
#SBATCH --error=ase_ose_logs_n2000/ase_ose_n2000_rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# ASE vs OSE VERTEX-WISE COVERAGE: n=2000 (no covariates)
# ============================================================================
# This is a pure network model (no covariates)
# Compares ASE and OSE vertex-wise coverage
#
# Expected time per rep: ~60-120 minutes (4x longer than n=1000)
# Resources: 1 CPU (serial), 32GB RAM
# ============================================================================

# Load R module
module load r/4.3.1

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scenario: ASE vs OSE (n=2000, no covariates)"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "========================================"

# Create log directory if needed
cd $SLURM_SUBMIT_DIR
mkdir -p ase_ose_logs_n2000

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as replication number
Rscript ase_ose_vertex_wise_n2000_100reps.R $SLURM_ARRAY_TASK_ID

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
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
fi
echo "========================================"

exit $EXIT_CODE
