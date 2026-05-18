# HPC Setup Instructions for n=2000 Parallel Simulation

## Overview
This guide explains how to update the `cgrdpg` package on HPC and run the n=2000 vertex-wise coverage simulation with parallel optimization.

## Step 1: Update the cgrdpg Package on HPC

### 1.1 SSH into HPC
```bash
ssh mle6@bigred200.uits.iu.edu
```

### 1.2 Navigate to cgrdpg directory
```bash
cd ~/cgrdpg
```

### 1.3 Pull latest changes from GitHub
```bash
# Make sure you're on the correct branch
git checkout claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN

# Pull latest changes (includes parallel optimization)
git pull origin claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN
```

### 1.4 Check that new files are present
```bash
# Verify parallel optimization files
ls -lh R/fisher_parallel.R
ls -lh simulations/vertex_wise_coverage_n2000_parallel_100reps.R
ls -lh simulations/submit_vertex_wise_n2000_parallel_100reps.slurm
```

### 1.5 Load R module
```bash
module load r/4.3.1
```

### 1.6 Install required dependencies
```bash
# Install foreach and doParallel packages (one-time setup)
R --vanilla << 'EOF'
install.packages(c("foreach", "doParallel"),
                 repos = "https://cloud.r-project.org",
                 lib = "~/R/x86_64-pc-linux-gnu-library/4.3")
EOF
```

### 1.7 Rebuild and install cgrdpg package
```bash
# From the cgrdpg directory
R CMD INSTALL --no-multiarch --with-keep.source .
```

Expected output should end with:
```
* DONE (cgrdpg)
```

### 1.8 Verify installation
```bash
R --vanilla << 'EOF'
library(cgrdpg)
# Check if parallel functions are available
if (exists("fit_grdpg_cov_parallel")) {
  cat("✓ Parallel functions loaded successfully!\n")
} else {
  cat("✗ ERROR: Parallel functions not found\n")
}
EOF
```

## Step 2: Prepare for Simulation

### 2.1 Navigate to simulations directory
```bash
cd ~/cgrdpg/simulations
```

### 2.2 Create log directory
```bash
mkdir -p vertex_wise_logs_n2000_parallel
```

### 2.3 Verify simulation files
```bash
ls -lh vertex_wise_coverage_n2000_parallel_100reps.R
ls -lh submit_vertex_wise_n2000_parallel_100reps.slurm
ls -lh aggregate_vertex_wise_n2000_results.R
```

## Step 3: Submit Array Job

### 3.1 Submit the job
```bash
sbatch submit_vertex_wise_n2000_parallel_100reps.slurm
```

You should see output like:
```
Submitted batch job XXXXXX
```

### 3.2 Check job status
```bash
# View your running/pending jobs
squeue -u mle6

# Count running jobs
squeue -u mle6 | grep vtx_n2000 | wc -l

# Check specific job
scontrol show job XXXXXX
```

### 3.3 Monitor progress
```bash
# Count completed result files
ls vertex_wise_results_n2000_parallel/*.rds 2>/dev/null | wc -l

# Watch progress (updates every 30 seconds)
watch -n 30 'echo "Completed: $(ls vertex_wise_results_n2000_parallel/*.rds 2>/dev/null | wc -l)/100"; squeue -u mle6 | grep vtx_n2000 | wc -l | xargs echo "Running/Pending:"'

# Press Ctrl+C to exit watch
```

### 3.4 Check for errors
```bash
# List non-empty error files (indicates problems)
find vertex_wise_logs_n2000_parallel/ -name "*.err" -type f ! -size 0

# Check specific error file
cat vertex_wise_logs_n2000_parallel/vertex_wise_n2000_rep001_*.err

# Check output file to see progress
tail -50 vertex_wise_logs_n2000_parallel/vertex_wise_n2000_rep001_*.out
```

## Step 4: After All Jobs Complete

### 4.1 Verify all 100 replications finished
```bash
ls vertex_wise_results_n2000_parallel/*.rds | wc -l
# Should output: 100
```

### 4.2 Run aggregation
```bash
module load r/4.3.1
Rscript aggregate_vertex_wise_n2000_results.R
```

### 4.3 Check aggregated results
```bash
ls -lh vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated_n2000.rds
ls -lh vertex_wise_results_n2000_parallel/vertex_wise_coverage_results_n2000.csv
```

## Step 5: Transfer Results to Local Machine

From your **local machine** (Mac), run:

```bash
# Transfer aggregated results
scp mle6@bigred200.uits.iu.edu:~/cgrdpg/simulations/vertex_wise_results_n2000_parallel/vertex_wise_coverage_aggregated_n2000.rds ~/Documents/GitHub/cgrdpg/simulations/

scp mle6@bigred200.uits.iu.edu:~/cgrdpg/simulations/vertex_wise_results_n2000_parallel/vertex_wise_coverage_results_n2000.csv ~/Documents/GitHub/cgrdpg/simulations/

# Or transfer entire results directory
scp -r mle6@bigred200.uits.iu.edu:~/cgrdpg/simulations/vertex_wise_results_n2000_parallel ~/Documents/GitHub/cgrdpg/simulations/
```

## Resource Allocation

The SLURM script requests:
- **CPUs**: 32 cores per job (for parallel model fitting)
- **Memory**: 128GB RAM (n=2000 requires more memory)
- **Time**: 96 hours (4 days) per job
- **Jobs**: 100 array jobs (one per replication)

### Expected Performance
- **Serial (n=2000)**: ~120 minutes per replication
- **Parallel (32 cores)**: ~25-35 minutes per replication
- **Speedup**: ~3-5x faster
- **Total time for 100 reps**: ~35-60 hours (1.5-2.5 days)

## Troubleshooting

### Issue: Package installation fails
```bash
# Check if R library directory exists
ls -ld ~/R/x86_64-pc-linux-gnu-library/4.3

# If not, create it
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.3

# Try installation again
R CMD INSTALL --no-multiarch --with-keep.source .
```

### Issue: Parallel functions not found
```bash
# Check if fisher_parallel.R was pulled
ls -lh R/fisher_parallel.R

# If missing, ensure you're on correct branch
git status
git checkout claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN
git pull
```

### Issue: Jobs failing with memory errors
```bash
# Check error logs
grep -i "memory\|oom\|killed" vertex_wise_logs_n2000_parallel/*.err

# If memory issues, increase memory in SLURM script:
# Change: #SBATCH --mem=128G
# To:     #SBATCH --mem=256G
```

### Issue: Jobs timing out
```bash
# Check timing from completed jobs
grep "Total runtime" vertex_wise_logs_n2000_parallel/*.out | tail -5

# If approaching 96 hours, increase time limit:
# Change: #SBATCH --time=96:00:00
# To:     #SBATCH --time=120:00:00
```

## Quick Reference Commands

```bash
# Submit job
sbatch submit_vertex_wise_n2000_parallel_100reps.slurm

# Check status
squeue -u mle6 | grep vtx_n2000

# Count completed
ls vertex_wise_results_n2000_parallel/*.rds | wc -l

# Cancel all jobs
scancel -u mle6 --name=vtx_n2000

# Cancel specific job
scancel JOBID

# View job details
scontrol show job JOBID

# Check node resources
sinfo -p general
```

## Notes

1. **Parallel Efficiency**: With 32 cores, expect ~4-5x speedup (vs theoretical 32x) due to:
   - Communication overhead between processes
   - Memory bandwidth limitations
   - Synchronization at each iteration

2. **Convergence**: Parallel version uses Jacobi-style updates which may differ slightly from serial Gauss-Seidel updates. Results should be within 1-2% SSE.

3. **Package Updates**: Whenever you pull new changes that modify R code, you must rebuild the package with `R CMD INSTALL`.

4. **Monitoring**: Check both `.out` and `.err` files regularly to catch issues early.
