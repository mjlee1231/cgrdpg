# 3D GRDPG Simulations

Vertex-wise coverage simulations for 3D Generalized Random Dot Product Graph (GRDPG) with signature matrix `S = diag(1, 1, -1)`.

## Overview

- **Latent space dimension**: d = 3
- **Signature matrix**: S = diag(1, 1, -1) (p=2, q=1)
- **Latent curve**: 3D helix
  - X₁ = 0.175 × cos(θ) + 0.530
  - X₂ = 0.175 × sin(θ) + 0.530
  - X₃ = 0.040 × θ/π + 0.315
  - θ ∈ [0, π] for i = 1,...,n
- **Edge probability range**: [0.281, 0.749] (width = 0.468)
- **Sample sizes**: n = 500, n = 1000
- **Replications**: 100 per sample size

## Parameters

| Parameter | n=500 | n=1000 | Description |
|-----------|-------|--------|-------------|
| maxit     | 10    | 10     | Maximum Fisher-scoring iterations |
| tol       | 0.05  | 0.05   | Convergence tolerance |
| tau       | 0.01  | 0.01   | Smoothing parameter for ψ function |
| p_cov     | 250   | 500    | Number of covariates (n/2) |
| ncores    | 16    | 16     | CPU cores for parallel fitting |

## Files

### Design and Testing
- `design_3d_curve_wider.R` - Optimize 3D helix curve for maximum edge prob range
- `3d_curve_design.rds` - Saved optimal curve design

### Simulation Scripts
- `vertex_wise_3d_n500_single.R` - Single replication for n=500
- `vertex_wise_3d_n1000_single.R` - Single replication for n=1000

### SLURM Submission
- `submit_3d_n500.slurm` - Submit 100 jobs for n=500 (24h, 64GB)
- `submit_3d_n1000.slurm` - Submit 100 jobs for n=1000 (48h, 128GB)

### Aggregation
- `aggregate_results_n500.R` - Aggregate and visualize n=500 results
- `aggregate_results_n1000.R` - Aggregate and visualize n=1000 results

## Directory Structure on HPC

After running SLURM jobs:

```
3d_simulations/
├── vertex_wise_3d_n500_single.R
├── vertex_wise_3d_n1000_single.R
├── submit_3d_n500.slurm
├── submit_3d_n1000.slurm
├── 3d_simulation_n500/          # Created automatically by SLURM
│   ├── logs/
│   │   ├── rep_1.out
│   │   ├── rep_1.err
│   │   ├── rep_2.out
│   │   └── ...
│   └── results_3d_n500/
│       ├── rep_001.rds
│       ├── rep_002.rds
│       └── ...
└── 3d_simulation_n1000/         # Created automatically by SLURM
    ├── logs/
    └── results_3d_n1000/
```

## Running Simulations

### On HPC

1. **Upload files to HPC:**
```bash
# From your local machine
scp vertex_wise_3d_n500_single.R vertex_wise_3d_n1000_single.R \\
    submit_3d_n500.slurm submit_3d_n1000.slurm \\
    user@hpc:/path/to/simulations/3d_simulations/
```

2. **Submit jobs:**
```bash
# On HPC
cd /path/to/simulations/3d_simulations/

# For n=500 (100 replications)
sbatch submit_3d_n500.slurm

# For n=1000 (100 replications)
sbatch submit_3d_n1000.slurm
```

3. **Monitor jobs:**
```bash
# Check job status
squeue -u $USER

# Check specific job array
squeue -j <job_id>

# View running/completed jobs
sacct -u $USER --format=JobID,JobName,State,Elapsed,MaxRSS
```

4. **Download results:**
```bash
# After jobs complete, download from HPC to local
scp -r user@hpc:/path/to/simulations/3d_simulations/3d_simulation_n500 ./
scp -r user@hpc:/path/to/simulations/3d_simulations/3d_simulation_n1000 ./
```

### Testing Locally

Before submitting to HPC, test with a single replication:

```bash
# Test n=500
Rscript vertex_wise_3d_n500_single.R 1

# Test n=1000
Rscript vertex_wise_3d_n1000_single.R 1
```

## Aggregating Results

After downloading results from HPC:

```bash
# Move results to expected location (if needed)
mv 3d_simulation_n500/results_3d_n500 ./
mv 3d_simulation_n1000/results_3d_n1000 ./

# Aggregate n=500 results
Rscript aggregate_results_n500.R

# Aggregate n=1000 results
Rscript aggregate_results_n1000.R
```

## Output Files

### Per Replication (*.rds in results_3d_n*)
Each replication saves:
- `in_ellipse_true`: Vertex-wise coverage indicators (TRUE G_in)
- `in_ellipse_plugin`: Vertex-wise coverage indicators (PLUGIN G_in)
- `SSE`: Sum of squared errors (alignment quality)
- `coverage_true/plugin`: Overall coverage rates for this replication
- `fit_time`, `coverage_time`: Timing information
- `converged`, `iterations`: Convergence diagnostics
- `n`, `p_cov`, `d`, `tau`, `maxit`, `tol`: Parameters used

### Aggregated Results
- `aggregated_results_3d_n*.rds` - All aggregated statistics across 100 reps
- `summary_3d_n*.csv` - Summary table (mean, SD, range)
- `coverage_boxplot_3d_n*.pdf` - Boxplot comparing TRUE vs PLUGIN coverage
- `vertex_coverage_3d_n*.pdf` - Vertex-wise coverage rates across curve

## Expected Runtime

Based on 2D simulations with similar settings:

- **n=500**: ~1-3 hours per replication
  - Model fitting: ~30-90 minutes (with 16 cores)
  - Coverage computation: ~10-30 minutes
  - **Time limit**: 24 hours (conservative)

- **n=1000**: ~4-12 hours per replication
  - Model fitting: ~2-6 hours (with 16 cores)
  - Coverage computation: ~1-2 hours
  - **Time limit**: 48 hours (conservative)

With parallelization (`fit_grdpg_cov_parallel`), these times should be on the lower end.

## Comparison with 1D and 2D

| Dimension | Signature Matrix | Latent Curve | Edge Prob Range | Width |
|-----------|-----------------|--------------|-----------------|-------|
| 1D        | diag(1)         | Sine wave    | varies          | ~0.84 |
| 2D        | diag(1,1)       | Circle/ellipse| [0.20, 0.66]   | 0.46  |
| **3D**    | **diag(1,1,-1)**| **Helix**    | **[0.28, 0.75]**| **0.47** |

The 3D GRDPG case provides:
- **Wider range than previous 2D**: Better for coverage estimation precision
- **Safer from boundaries**: min=0.28 and max=0.75 are far from 0 and 1
- **Richer geometry**: Helix structure in 3D with indefinite signature matrix
- **Realistic setting**: Tests GRDPG with both positive and negative eigenvalues

## Troubleshooting

### Out of Memory
If jobs fail with OOM errors:
- Increase `--mem` in SLURM script (currently 64GB for n=500, 128GB for n=1000)
- Reduce `ncores` if memory per core is limiting

### Timeout
If jobs hit time limit:
- Check logs to see where it's spending time
- May need to reduce `maxit` or increase time limit
- Consider using only parallel fitting version

### Convergence Issues
If many replications don't converge:
- Check `tau` value (currently 0.01)
- May need to adjust `tol` or `maxit`
- Examine failed replications' logs

### Missing Parallel Function
If `fit_grdpg_cov_parallel` not available:
- Scripts will automatically fall back to serial `fit_grdpg_cov`
- Will be slower but should still work
- Make sure cgrdpg package is up to date

## Contact

For questions about these simulations, refer to the main simulation documentation or check logs for specific error messages.
