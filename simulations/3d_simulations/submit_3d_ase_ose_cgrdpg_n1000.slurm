#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=3d_ase_ose_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=1-100
#SBATCH --output=logs/rep%03a_%j.out
#SBATCH --error=logs/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# 3D GRDPG: ASE vs OSE vs cgrdpg vertex-wise coverage
# n=1000, 100 replications (1 rep per node), p_cov=500, S=diag(1,1,-1)

module load r/4.3.1

mkdir -p logs
mkdir -p results_3d_ase_ose_n1000

echo "========================================"
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Methods: cgrdpg-TRUE/PLUGIN, ASE-TRUE/PLUGIN, OSE-TRUE/PLUGIN"
echo "Replication: $SLURM_ARRAY_TASK_ID/100"
echo "Resources: 32 CPUs, 128GB RAM"
echo "========================================"

cd $SLURM_SUBMIT_DIR

START_TIME=$(date +%s)

Rscript ase_ose_vertex_wise_3d_n1000.R $SLURM_ARRAY_TASK_ID

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
