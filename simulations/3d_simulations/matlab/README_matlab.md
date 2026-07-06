# MATLAB Implementation of CGRDPG with fminunc

This directory contains MATLAB implementations of CGRDPG (Covariate-assisted Generalized Random Dot Product Graph) optimization using MATLAB's `fminunc` function.

## Motivation

The R implementation using Fisher scoring showed numerical instability in 3D simulations (n=1000, 2000):
- **OSE**: Extreme outliers (SSE up to 1.5 million)
- **CGRDPG (Fisher)**: Occasional convergence failures

MATLAB's `fminunc` offers:
1. **Multiple robust algorithms**: quasi-Newton (BFGS), trust-region
2. **Automatic gradient verification**: Checks analytical gradients against numerical ones
3. **Better line search**: Strong Wolfe conditions for step size
4. **Superior numerical stability**: Mature optimization toolkit
5. **Detailed diagnostics**: First-order optimality, condition numbers, etc.

## Directory Structure

```
matlab/
├── core/                                    # Core implementation
│   ├── fit_grdpg_fminunc_surrogate.m       # Main optimization function
│   └── psi_functions.m                      # Smoothed log-likelihood functions
├── test_n1000_p500.m                        # Single test (n=1000, p_cov=500)
├── test_gradient.m                          # Gradient correctness checker
├── quick_test_surrogate.m                   # Quick test (n=500, p_cov=250)
├── submit_test_n1000_p500.slurm            # HPC submission script
└── [other test and legacy files]
```

## Core Files (`core/` folder)

### `fit_grdpg_fminunc_surrogate.m`
Main optimization function using surrogate log-likelihood (matches R cgrdpg package):

**Features:**
- **Initialization**: ASE (Adjacency Spectral Embedding)
- **Objective**: Surrogate log-likelihood with psi functions
- **Optimizer**: MATLAB's fminunc with quasi-Newton
- **Outputs**: Optimized latent positions X and covariate coefficients Z

**Usage:**
```matlab
addpath('core');
[X_opt, Z_opt, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);
```

### `psi_functions.m`
Smoothed log-likelihood helper functions:

- `psi(x)`: Smoothed log-odds ≈ log(x/(1-x))
- `Psi(x)`: Integral of psi(x)
- `dpsi(x)`: Derivative of psi(x)

Uses quadratic patches for smooth extension outside [0,1].

## Test Files

### `test_n1000_p500.m`
Single replication test matching R specification:
- n=1000 nodes, p_cov=500 covariates, d=3 dimensions
- Parametric latent curve: X = [sin(2πt)+0.6, cos(2πt)+0.6, cos(4πt)]
- Compares ASE vs Surrogate optimization
- Verifies objective function components

### `test_gradient.m`
Gradient correctness verification using finite differences:
- Tests analytical gradient against numerical approximation
- Reports max relative error
- Critical for debugging optimization issues

### `quick_test_surrogate.m`
Fast test for development (n=500, p_cov=250)

## Quick Start

### 1. Run 100 replications (recommended)

```matlab
% In MATLAB:
cd /Users/jin/Documents/GitHub/cgrdpg/simulations/3d_simulations/matlab
run_100_replications
```

This will:
- Run 100 independent replications
- Match R simulation setup (n=1000, d=3)
- Save summary statistics and individual results
- Take approximately 30-60 minutes (depends on CPU)

Output files:
- `results/fminunc_100reps_n1000_summary.mat` - Summary statistics
- `results/fminunc_100reps_n1000_detailed.mat` - All replications
- `results/replications/rep_*.mat` - Individual replication files

### 2. Run a single test (quick check)

```matlab
test_3d_simulation
```

This will:
- Generate synthetic data (n=1000, p_cov=1000, d=3)
- Run single fminunc optimization
- Compare with ASE baseline
- Save results to `matlab_3d_test_results.mat`

### 3. Compare with R results

```matlab
compare_matlab_vs_r
```

This creates distribution plots and statistical comparisons between MATLAB and R.

### 4. Configure fminunc options

```matlab
% Quasi-Newton (default, fast)
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'iter', ...
    'MaxIterations', 100);

% Trust-region (more stable, requires gradient)
options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'bfgs', ...
    'Display', 'iter');
```

### 3. Compare with R results

First, export R results to MATLAB format:

```r
# In R:
library(R.matlab)
results <- readRDS("simulations/3d_simulations/results_3d_ase_ose_cgrdpg_n1000/rep_001.rds")
writeMat("matlab/r_results_rep001.mat", results=results)
```

Then in MATLAB:

```matlab
compare_with_r
```

## Function Usage

### Basic usage

```matlab
% Add core folder to path
addpath('core');

% Load your data
A = ...; % n x n adjacency matrix
B = ...; % p_cov x n covariate matrix
d = 3;   % embedding dimension
p = 2;   % positive signatures (q = d - p)
tau = 0.001;

% Fit model
[X_opt, Z_opt, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau);

% Check convergence
if exitflag > 0
    fprintf('Converged successfully!\n');
    fprintf('Final objective: %.6f\n', fval);
    fprintf('Iterations: %d\n', output.iterations);
else
    fprintf('Warning: Did not converge (exitflag = %d)\n', exitflag);
end
```

### Advanced options

```matlab
% Custom options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'iter-detailed', ...  % Verbose output
    'MaxIterations', 200, ...
    'OptimalityTolerance', 1e-8, ... % Stricter tolerance
    'CheckGradients', true);         % Verify gradients

[X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc(A, B, d, p, tau, options);
```

## Understanding Exit Flags

- `exitflag = 1`: First-order optimality satisfied (converged)
- `exitflag = 2`: Step size below threshold (likely converged)
- `exitflag = 0`: Maximum iterations reached
- `exitflag = -1`: Stopped by output function
- `exitflag = -3`: Objective unbounded below (problem with data/model)

## Expected Performance

Based on preliminary tests:

| Method | Typical SSE | Time | Stability |
|--------|-------------|------|-----------|
| ASE | ~70-80 | <1s | Very stable |
| OSE (R) | ~130-500k | <1s | **Unstable** |
| Fisher (R) | ~5-250 | 30-60s | Occasional failures |
| fminunc | ~5-20 | 10-30s | **Expected: More stable** |

## Troubleshooting

### Gradient verification failed
```matlab
options.CheckGradients = true;
% Run and check the gradient error - should be < 1e-6
```

### Slow convergence
```matlab
% Try trust-region algorithm
options.Algorithm = 'trust-region';
% Or increase tolerance
options.OptimalityTolerance = 1e-5;
```

### Numerical issues (NaN/Inf)
```matlab
% Increase tau for more stability
tau = 0.01;  % instead of 0.001
% Check condition of Hessian
options.Display = 'iter-detailed';
```

## Running on HPC (Recommended!)

### Why HPC?

Running on HPC is **much faster** than local execution:
- **Parallel processing**: 48 cores run replications simultaneously
- **High-performance CPUs**: Faster per-iteration computation
- **Uninterrupted**: Background execution for hours
- **Expected speedup**: 30-60 min → ~1-2 hours for 100 reps

### Setup on HPC

1. **Check MATLAB availability:**

```bash
cd ~/cgrdpg/simulations/3d_simulations/matlab
bash check_hpc_matlab.sh
```

2. **Submit the job(s):**

```bash
# Make sure you're in the matlab directory
cd ~/cgrdpg/simulations/3d_simulations/matlab

# Create log/error directories
mkdir -p logs errors results

# Option 1: Submit both n=1000 and n=2000 (recommended)
bash submit_both.sh

# Option 2: Submit individually
sbatch submit_matlab_100reps.slurm           # n=1000
sbatch submit_matlab_100reps_n2000.slurm    # n=2000
```

3. **Monitor progress:**

```bash
# Check job status
squeue -u mle6

# Watch the log file in real-time
tail -f logs/matlab_100reps_*.out

# Check progress
ls -lh results/replications_parallel/ | wc -l  # Count completed reps
```

4. **After completion:**

```bash
# Check results
ls -lh results/

# Download results to local machine (from your Mac)
cd ~/Documents/GitHub/cgrdpg/simulations/3d_simulations/matlab
rsync -avz mle6@bigred200.uits.iu.edu:~/cgrdpg/simulations/3d_simulations/matlab/results/ ./results/
```

### HPC Job Configuration

**For n=1000** (`submit_matlab_100reps.slurm`):
- **CPUs**: 48 cores (parallel parfor execution)
- **Memory**: 192GB
- **Time limit**: 4 hours
- **Expected runtime**: 1-2 hours

**For n=2000** (`submit_matlab_100reps_n2000.slurm`):
- **CPUs**: 48 cores (parallel parfor execution)
- **Memory**: 256GB (increased for larger network)
- **Time limit**: 8 hours
- **Expected runtime**: 2-4 hours

### Troubleshooting HPC

**MATLAB module not found:**
```bash
module spider matlab
module load matlab/R2023a  # Use available version
```

**Out of memory:**
```bash
# Edit submit_matlab_100reps.slurm
#SBATCH --mem=256G  # Increase from 192G
```

**Job timed out:**
```bash
# Edit submit_matlab_100reps.slurm
#SBATCH --time=08:00:00  # Increase from 4 hours
```

## Next Steps

1. ✅ **Run 100 reps on HPC** (using instructions above)
2. **Download and compare results** with R
3. **Generate SSE histograms** using `compare_matlab_vs_r.m`
4. **Test on real data**: LastFM dataset (n=1699, d=1 and d=4)

## References

- MATLAB fminunc documentation: https://www.mathworks.com/help/optim/ug/fminunc.html
- MATLAB Parallel Computing: https://www.mathworks.com/help/parallel-computing/
- Quasi-Newton methods: Nocedal & Wright, "Numerical Optimization"
- GRDPG: Athreya et al. (2017), "Statistical inference on random dot product graphs"
