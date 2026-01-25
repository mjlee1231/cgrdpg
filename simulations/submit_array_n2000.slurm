#!/bin/bash
#SBATCH --account=r00357
#SBATCH --job-name=grdpg_n2000_array
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-100
#SBATCH --output=results_n2000/output_n2000_rep%a.out
#SBATCH --error=results_n2000/output_n2000_rep%a.err
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
echo "Scenario: n=2000, p_cov=1000"
echo "========================================"
echo

# Navigate to submission directory and create results directory
cd $SLURM_SUBMIT_DIR
mkdir -p results_n2000

# Run the R script (SLURM_ARRAY_TASK_ID is automatically set)
Rscript single_rep_n2000_ellipse.R

echo
echo "========================================"
echo "Job finished at: $(date)"
echo "========================================"
