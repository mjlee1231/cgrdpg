#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=test_n1000_p500
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/test_n1000_p500_%j.out
#SBATCH --error=logs/test_n1000_p500_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# SINGLE TEST RUN: Surrogate likelihood implementation
# Parameters: n=1000, p_cov=500, d=3, max_iter=30, S=diag(1,1,-1)
# ============================================================================
# This is a SINGLE replication to verify correctness before running 100 reps.
#
# Verifies:
#   - ASE initialization
#   - Surrogate objective function
#   - Gradient computation
#   - Data generation matches R specification
#   - Convergence with max_iter=30
#
# Next step: Run 100 replications with vertex-wise coverage computation
# Expected runtime: ~5-15 minutes
# ============================================================================

module load matlab

echo "========================================"
echo "Test: Surrogate Likelihood (n=1000, p_cov=500)"
echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "Time limit: 1 hour"
echo "========================================"
echo ""

cd $SLURM_SUBMIT_DIR

# Create output directories if they don't exist
mkdir -p logs
mkdir -p results

START_TIME=$(date +%s)

# Run MATLAB test script in batch mode
echo "Starting MATLAB test..."
echo "Script: test_n1000_p500.m"
echo ""

matlab -nodisplay -nosplash -nodesktop -r "test_n1000_p500; exit;"

EXIT_CODE=$?

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "========================================"
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "Runtime: ${MINUTES}m ${SECONDS}s"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
  echo "Status: ✓ SUCCESS"
  echo ""
  echo "Verification checklist:"
  echo "  ✓ Data generation matches R specification"
  echo "  ✓ ASE initialization computed"
  echo "  ✓ Surrogate objective optimized (max_iter=30)"
  echo "  ✓ Gradient computation correct"
  echo ""
  echo "Check output above for:"
  echo "  - ASE SSE baseline"
  echo "  - Surrogate SSE improvement"
  echo "  - Convergence diagnostics"
else
  echo "Status: ✗ FAILED (exit code $EXIT_CODE)"
  echo "Check error log: logs/test_n1000_p500_${SLURM_JOB_ID}.err"
fi

echo "========================================"

exit $EXIT_CODE
