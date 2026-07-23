#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=matlab_2d_no_info_n1000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-100
#SBATCH --output=logs_matlab_2d_no_info_n1000/rep%03a_%j.out
#SBATCH --error=logs_matlab_2d_no_info_n1000/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# 2D RDPG coverage simulation with NO COVARIATE INFORMATION
# n=1000, p_cov=500, Z0=0, B=noise only
# Methods: cgrdpg, ASE, OSE

# Create output directories
mkdir -p logs_matlab_2d_no_info_n1000
mkdir -p results_2d_no_info_n1000

# Load MATLAB module
module load matlab/2025a

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "NO COVARIATE INFO: Z0=0, B=noise"
echo "Sample size: n=1000, p_cov=500"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

# Run single replication
matlab -batch "coverage_2d_no_info_single_rep_n1000(${SLURM_ARRAY_TASK_ID})"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
