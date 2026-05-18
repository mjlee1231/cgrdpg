#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=vtx_n3000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --array=1-100
#SBATCH --output=vertex_wise_logs_n3000_rcpp/vertex_wise_n3000_rep%03a_%j.out
#SBATCH --error=vertex_wise_logs_n3000_rcpp/vertex_wise_n3000_rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# VERTEX-WISE COVERAGE: n=3000 with PARALLEL + RCPP OPTIMIZATION
# ============================================================================
# Expected performance:
#   - Parallel optimization: ~4x speedup (64 cores)
#   - Rcpp optimization: ~2x additional speedup
#   - Total: ~6-8x faster than serial
#   - Estimated time per rep: 40-50 minutes (vs ~3 hours serial)
#   - Total for 100 reps: ~60-80 hours (2.5-3.5 days)
# ============================================================================

# Load R module
module load r/4.3.1

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scenario: Vertex-wise Coverage - PARALLEL + RCPP (n=3000, p_cov=1500)"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: $SLURM_CPUS_PER_TASK CPUs, $SLURM_MEM_PER_NODE MB RAM"
echo "Optimization: Parallel (4x) + Rcpp C++ (2x) = 6-8x total speedup"
echo "========================================"

# Create log directory if needed
cd $SLURM_SUBMIT_DIR
mkdir -p vertex_wise_logs_n3000_rcpp

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as replication number
Rscript vertex_wise_coverage_n3000_rcpp_100reps.R $SLURM_ARRAY_TASK_ID

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
