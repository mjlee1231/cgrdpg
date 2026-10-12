# MATLAB Dimension Selection for LastFM Data

MATLAB implementation of dimension selection using `fminunc` with trust-region algorithm.

## Why MATLAB?

**Expected advantages over R version:**
- **Faster optimization**: MATLAB's `fminunc` with exact Hessian is highly optimized
- **Better convergence**: Trust-region method with exact derivatives
- **Efficiency**: 2-3x speedup expected (R took 3.96 hours for d=1)

## Files

- `dimension_selection_single.m` - Main script for single dimension
- `submit_matlab_dimension_selection.slurm` - HPC SLURM submission script
- `prepare_data_for_matlab.R` - Convert R data to MATLAB format
- `core/` - Core MATLAB functions (copied from simulations)

## Setup

### 1. Prepare Data (One-time)

```bash
# On local machine or HPC
cd /path/to/LastFM/matlab_dimension_selection
Rscript prepare_data_for_matlab.R
```

This creates `CovariatesMatrix.mat` from the R data file.

### 2. Copy CSV Files

```bash
# Create symbolic links to avoid duplication
ln -s ../lastfm_asia_edges.csv .
ln -s ../lastfm_asia_target.csv .
```

## Running on HPC

### Submit All Dimensions (d=1,2,3,4,5)

```bash
# SSH to HPC
ssh mle6@bigred200.uits.iu.edu

# Navigate to directory
cd ~/cgrdpg/LastFM/matlab_dimension_selection

# Submit array job
module load matlab
sbatch submit_matlab_dimension_selection.slurm

# Check status
squeue -u mle6

# Monitor logs
tail -f logs_matlab/dimension_d1_*.out
```

### Results

Results are saved to `results_matlab_maxit30_tol0.1/results_d{1-5}.mat`

Each `.mat` file contains:
- `d` - dimension
- `neg_log_lik` - negative log pseudo-likelihood
- `n_params` - number of parameters
- `BIC` - Bayesian Information Criterion
- `AIC` - Akaike Information Criterion
- `converged` - convergence status
- `iters` - number of iterations
- `time_sec` - elapsed time in seconds
- `latent_positions` - fitted X matrix
- `covariate_effects` - fitted Z matrix

## Expected Performance

Based on R version timing (d=1 took 3.96 hours with sequential processing):

| Dimension | R (sequential) | MATLAB (expected) |
|-----------|---------------|-------------------|
| d=1       | 3.96 hours    | **1-2 hours**     |
| d=2       | ~4-5 hours    | **1.5-2.5 hours** |
| d=3       | ~5-6 hours    | **2-3 hours**     |
| d=4       | ~6-7 hours    | **2.5-3.5 hours** |
| d=5       | ~7-8 hours    | **3-4 hours**     |

**Total time**: All 5 dimensions complete in **8 hours** (vs 20+ hours for R)

## Compare with R Results

```matlab
% Load MATLAB results
load('results_matlab_maxit30_tol0.1/results_d1.mat');
fprintf('MATLAB d=1: BIC=%.2f, Time=%.2f hours\n', result.BIC, result.time_sec/3600);

% Compare with R results (if available)
% R d=1: BIC=224529.54, Time=3.96 hours (maxit=30, tol=0.1)
```

## Troubleshooting

If MATLAB module not found:
```bash
module spider matlab
module load matlab/<version>
```

If memory issues:
- Increase `--mem` in SLURM script (currently 64G)
- Use `--mem=128G` for larger dimensions
