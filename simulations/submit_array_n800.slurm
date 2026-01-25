#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=grdpg_n800_array
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-3
#SBATCH --output=results_n800/output_n800_rep%a.out
#SBATCH --error=results_n800/output_n800_rep%a.err
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
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Scenario: n=800, p_cov=400"
echo "========================================"
echo

# Navigate to submission directory and create results directory
cd $SLURM_SUBMIT_DIR
mkdir -p results_n800

# Run the R script (SLURM_ARRAY_TASK_ID is automatically set)
Rscript single_rep_n800_ellipse.R

echo
echo "========================================"
echo "Job finished at: $(date)"
echo "========================================"
