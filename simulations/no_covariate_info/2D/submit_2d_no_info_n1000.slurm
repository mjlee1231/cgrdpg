#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=2d_noinfo_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-100
#SBATCH --output=logs_2d_no_info_n1000/rep%03a_%j.out
#SBATCH --error=errors_2d_no_info_n1000/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# 2D No Covariate Info: cgrdpg vs ASE vs OSE vertex-wise coverage
# Z0 = 0, B = pure noise (no signal)
# n=1000, 100 replications, p_cov=500, d=2 (RDPG)

module load r/4.3.1

# Create separate directories for results, logs, and errors
mkdir -p results_2d_no_info_n1000
mkdir -p logs_2d_no_info_n1000
mkdir -p errors_2d_no_info_n1000

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Scenario: Z0 = 0, B = pure noise"
echo "Methods: cgrdpg-TRUE/PLUGIN, ASE-TRUE/PLUGIN, OSE-TRUE/PLUGIN"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: 16 CPUs, 64GB RAM"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

Rscript vertex_wise_coverage_2d_no_info_n1000.R $SLURM_ARRAY_TASK_ID

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
