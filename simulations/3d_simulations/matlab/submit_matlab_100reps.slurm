#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_fminunc_100
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=04:00:00
#SBATCH --output=logs/matlab_100reps_n1000_%j.out
#SBATCH --error=errors/matlab_100reps_n1000_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# ============================================================================
# MATLAB fminunc: 100 replications (3D CGRDPG, n=1000)
# ============================================================================
# Uses 48 cores for parallel execution
# Expected time: 1-2 hours (much faster than serial ~30-60 min)
# ============================================================================

module load matlab

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory: 192GB"
echo "Time limit: 4 hours"
echo "========================================"

cd $SLURM_SUBMIT_DIR

# Create output directories if they don't exist
mkdir -p results/n1000/replications
mkdir -p logs
mkdir -p errors

START_TIME=$(date +%s)

# Run MATLAB script in batch mode
echo "Starting MATLAB..."
matlab -nodisplay -nosplash -nodesktop -r "run_100_replications_parallel; exit;"

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
  echo "100 replications completed"
  echo "Check results in: results/n1000/"
else
  echo "Status: FAILED (exit code $EXIT_CODE)"
  echo "Check error log: errors/matlab_100reps_n1000_${SLURM_JOB_ID}.err"
fi
echo "========================================"

exit $EXIT_CODE
