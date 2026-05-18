# GRDPG Simulation Instructions

## Overview

This simulation tests a rank-three generic GRDPG with signature (2,1), where latent positions are drawn from a parametric curve in R³.

### Latent Curve Parameterization

```
X(t) = [0.15*sin(2*pi*t) + 0.6,
        0.15*cos(2*pi*t) + 0.6,
        0.15*cos(4*pi*t)]^T
```

where `t ∈ [0,1]` and the n latent positions are obtained by `t = i/n` for `i ∈ [n]`.

### Simulation Scenarios

1. **Scenario 1**: n = 100, p_cov = 50
2. **Scenario 2**: n = 500, p_cov = 250

Each scenario is repeated 100 times.

### Model Parameters

- **Signature**: (p=2, q=1)
- **Embedding dimension**: d = 3
- **Covariate loading**: Z₀ ~ N(0, 0.5²)
- **Covariate noise**: W ~ N(0, 1)

## Files

- **`simulations_latent_curve.R`**: Main simulation script
- **`test_simulation.R`**: Quick test with 2 repetitions

## How to Run

### Prerequisites

Make sure the `cgrdpg` package is installed:

```r
# In R:
install.packages("devtools")
devtools::document("cgrdpg")
devtools::install("cgrdpg")
```

### Quick Test (Recommended First)

Run a quick test with just 2 repetitions to verify everything works:

```bash
Rscript test_simulation.R
```

This should complete in under a minute and verify that the code runs without errors.

### Full Simulation

Run the complete simulation (100 repetitions × 2 scenarios):

```bash
Rscript simulations_latent_curve.R
```

**Expected runtime**:
- Scenario 1 (n=100): ~10-30 seconds per repetition
- Scenario 2 (n=500): ~1-5 minutes per repetition
- Total: ~2-10 hours depending on your machine

You can also run it from within R:

```r
source("simulations_latent_curve.R")
```

## Output Files

The script generates:

1. **`simulation_results_latent_curve.rds`**: Full results including all SSE values and timing data
2. **`simulation_summary_latent_curve.csv`**: Summary table with mean/SD of SSE and mean runtime

## Output Format

The script reports:

| Scenario | N_Valid | Mean_SSE | SD_SSE | Mean_Runtime_sec |
|----------|---------|----------|--------|------------------|
| n=100, p_cov=50 | 100 | ... | ... | ... |
| n=500, p_cov=250 | 100 | ... | ... | ... |

Where:
- **Mean_SSE**: Average sum of squared errors after Procrustes alignment
- **SD_SSE**: Standard deviation of SSE across repetitions
- **Mean_Runtime_sec**: Average fitting time per repetition

## Implementation Details

### Data Generation

For each repetition:

1. Generate true latent positions `X₀` from the parametric curve
2. Compute probability matrix `P = X₀ S X₀ᵀ` where `S = diag(1, 1, -1)`
3. Clip probabilities to `[1e-3, 1-1e-3]` for numerical stability
4. Sample symmetric adjacency matrix `A ~ Bernoulli(P)` with no self-loops
5. Generate covariate loadings `Z₀ ~ N(0, 0.5²)` of size `p_cov × 3`
6. Generate covariate matrix `B = Z₀X₀ᵀ + W` where `W ~ N(0, 1)`

### Model Fitting

Fit the GRDPG model using:

```r
fit_grdpg_cov(A, B, d=3, p=2, q=1, maxit=5, tol=1e-2, tau=0.05)
```

### Error Computation

SSE is computed as:

```
SSE = ||X̂Q - X₀||²_F
```

where `Q` is the optimal orthogonal Procrustes transformation:
- Compute SVD: `UΣVᵀ = X̂ᵀX₀`
- Set `Q = UVᵀ`
- Align: `X_aligned = X̂Q`

This accounts for rotational/reflection ambiguity in the embedding space.

## Customization

You can modify simulation parameters by editing `simulations_latent_curve.R`:

- **Number of repetitions**: Change `n_reps = 100` in the `run_simulation()` calls
- **Model parameters**: Adjust `maxit`, `tol`, `tau` in the `run_simulation()` calls
- **Covariate noise levels**: Modify `sd` parameters in the `Z0` and `B` generation
- **Additional scenarios**: Add more `run_simulation()` calls with different `n` values

## Troubleshooting

If you encounter errors:

1. **Package not loaded**: Make sure `cgrdpg` is installed via `devtools::install("cgrdpg")`
2. **Memory issues**: For n=500, you may need at least 4-8GB RAM
3. **Convergence failures**: Adjust `maxit`, `tau`, or initialization parameters
4. **Slow runtime**: Reduce `n_reps` for testing, or run scenarios separately

## Results Interpretation

- **Lower SSE**: Better recovery of true latent positions
- **Lower SD**: More stable/consistent estimation across random instances
- **SSE typically increases with n**: More parameters to estimate, but relative error may decrease
