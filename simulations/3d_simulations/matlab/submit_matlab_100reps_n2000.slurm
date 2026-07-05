#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_fminunc_n2000
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=logs/matlab_100reps_n2000_%j.out
#SBATCH --error=errors/matlab_100reps_n2000_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# MATLAB fminunc: 100 replications (3D CGRDPG, n=2000)
# ============================================================================
# Uses 48 cores for parallel execution
# Expected time: 2-4 hours (larger network, more computation)
# Increased memory to 256GB for n=2000
# ============================================================================

module load matlab

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Network size: n = 2000"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory: 256GB"
echo "Time limit: 8 hours"
echo "========================================"

cd $SLURM_SUBMIT_DIR

# Create output directories if they don't exist
mkdir -p results/n2000/replications
mkdir -p logs
mkdir -p errors

START_TIME=$(date +%s)

# Run MATLAB script in batch mode
echo "Starting MATLAB for n=2000..."
matlab -nodisplay -nosplash -nodesktop -r "run_100_replications_n2000_parallel; exit;"

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
  echo "100 replications (n=2000) completed"
  echo "Check results in: results/n2000/"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
  echo "Check error log: errors/matlab_100reps_n2000_${SLURM_JOB_ID}.err"
fi
echo "========================================"

exit $EXIT_CODE
