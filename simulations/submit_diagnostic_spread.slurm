#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=diag_spread_%a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=248G
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/diag_spread_rep%a_%j.out
#SBATCH --error=logs/diag_spread_rep%a_%j.err
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
echo "Scenario: 2D SPREAD Oracle Alignment"
echo "Replication: $SLURM_ARRAY_TASK_ID/5"
echo "Resources: 128 CPUs, 248GB RAM, EXCLUSIVE node"
echo "CPUs available for R: $SLURM_CPUS_PER_TASK"
echo "========================================"

# Run the simulation with the replication number
cd $SLURM_SUBMIT_DIR

# Record start time
START_TIME=$(date +%s)

# Run R script with array task ID as replication number
Rscript diagnostic_2d_spread_oracle_5reps.R $SLURM_ARRAY_TASK_ID

# Record end time and compute duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo "========================================"
echo "Job finished at: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo "========================================"
