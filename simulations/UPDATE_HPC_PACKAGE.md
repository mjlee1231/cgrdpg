# Instructions: Update cgrdpg Package on HPC

## Step 1: SSH to HPC
```bash
ssh <your_username>@bigred200.uits.iu.edu
# or your HPC login node
```

## Step 2: Navigate to your cgrdpg directory
```bash
cd /path/to/your/cgrdpg
# Example: cd ~/cgrdpg or cd $HOME/cgrdpg
```

## Step 3: Pull latest changes
```bash
git fetch origin
git checkout claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN
git pull origin claude/latent-position-coverage-01CskqEusYY2gcrA7a6nEHSN
```

## Step 4: Load R module
```bash
module load r/4.3.1
```

## Step 5: Reinstall the package
```bash
R
```

Then in R:
```r
# Remove old version
remove.packages("cgrdpg")

# Install updated version from current directory
install.packages(".", repos = NULL, type = "source")

# Verify installation
library(cgrdpg)

# Check that fit_grdpg_cov exists
?fit_grdpg_cov

# Exit R
q()
```

## Step 6: Navigate to simulations directory and submit jobs
```bash
cd simulations

# Check that the new SLURM scripts exist
ls -lh submit_diagnostic_*.slurm

# Submit the array jobs
sbatch submit_diagnostic_tight.slurm
sbatch submit_diagnostic_spread.slurm
```

## Step 7: Monitor job status
```bash
# Check job queue
squeue -u $USER

# Check specific job details
scontrol show job <job_id>

# View output in real-time (once job starts)
tail -f diag_tight_rep1_*.out
tail -f diag_spread_rep1_*.out
```

## Expected Results
- Each scenario will launch 5 array jobs (one per replication)
- Total of 10 jobs running in parallel
- Each job outputs to: `diag_{tight|spread}_rep{1-5}_{jobid}.out`
- Results saved to: `diagnostic_2d_{tight|spread}_oracle_rep{1-5}.rds`
- Runtime per rep will be measured and reported

## Key Changes in This Update
1. **Adaptive step sizes** in `fit_grdpg_cov()`:
   - ls_beta = 0.8 for iterations 1-8
   - ls_beta = 0.4 for iterations 9+
2. **Array job support**: Each replication runs independently on separate node
3. **High-performance resources**: 32GB RAM, 16 CPUs per job
