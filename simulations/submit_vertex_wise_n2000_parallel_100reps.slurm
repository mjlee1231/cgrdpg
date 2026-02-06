#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=vtx_n2000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=96:00:00
#SBATCH --array=1-100
#SBATCH --output=vertex_wise_logs_n2000_parallel/vertex_wise_n2000_rep%03a_%j.out
#SBATCH --error=vertex_wise_logs_n2000_parallel/vertex_wise_n2000_rep%03a_%j.err
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
echo "Scenario: Vertex-wise Coverage - PARALLEL (n=2000, p_cov=1000)"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: $SLURM_CPUS_PER_TASK CPUs, $SLURM_MEM_PER_NODE MB RAM"
echo "Parallel optimization: ~3-6x speedup expected"
echo "========================================"

# Run the simulation with the replication number
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as replication number
Rscript vertex_wise_coverage_n2000_parallel_100reps.R $SLURM_ARRAY_TASK_ID

# Record end time and compute duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
