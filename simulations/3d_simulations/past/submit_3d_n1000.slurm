#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=3d_n1000_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=1-100
#SBATCH --output=3d_simulation_n1000/logs/rep%03a_%j.out
#SBATCH --error=3d_simulation_n1000/logs/rep%03a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mle6@iu.edu

# 3D GRDPG Simulation: n=1000, 100 replications
# Signature matrix: S = diag(1,1,-1)
# Each job runs one replication

module load r/4.3.1

# Create directory structure on HPC
mkdir -p 3d_simulation_n1000/logs
mkdir -p 3d_simulation_n1000/results_3d_n1000

# Change to simulation directory
cd 3d_simulation_n1000

# Run simulation for this array task ID
Rscript ../vertex_wise_3d_n1000_single.R ${SLURM_ARRAY_TASK_ID}
