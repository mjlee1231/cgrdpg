#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=vtx_opt_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-100
#SBATCH --output=vertex_wise_logs_n500_optimized/vertex_wise_spread_rep%03a_%j.out
#SBATCH --error=vertex_wise_logs_n500_optimized/vertex_wise_spread_rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# Load R module
module load r/4.3.1

# Print job information
echo "========================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scenario: Vertex-wise Coverage - 2D SPREAD (OPTIMIZED r=0.28, c=0.42)"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: 16 CPUs, 32GB RAM"
echo "========================================"

# Run the simulation with the replication number
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as replication number
Rscript vertex_wise_coverage_spread_100reps.R $SLURM_ARRAY_TASK_ID

# Record end time and compute duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
