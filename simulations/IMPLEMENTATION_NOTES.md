# Implementation Notes: Latent Position Coverage Testing

## Overview

This implementation tests the asymptotic normality and coverage probability of vertex-wise confidence sets for latent position vectors in the GRDPG model with covariates.

## Theoretical Background

### Model

The model assumes:
- Adjacency matrix: A ~ GRDPG(X, S) where S is the signature matrix
- Covariates: B = ZX^T + noise
- Latent positions follow a smooth curve in R^d

### Asymptotic Theory

For vertex i, the estimated position x̂_i satisfies:

```
√(n+p) * (x̂_i - x_{0i}) → N(0, Σ_i)
```

where the asymptotic covariance is:

```
Σ_i = G_in^{-1}
```

with

```
G_in = (1/(n+p)) * [Σ_{j=1}^n Ψ'(x_{0i}^T y_{0j}) y_{0j} y_{0j}^T + Σ_{l=1}^p z_{0l} z_{0l}^T]
```

Here:
- Ψ'(·) is the derivative of the smoothed log-likelihood (dpsi function)
- y_{0j} = S x_{0j} (signed latent positions)
- z_{0l} are the covariate loading vectors

### Confidence Sets

A (1-α)% confidence interval for coordinate k of vertex i is:

```
x̂_i[k] ± z_{α/2} * √(Σ_i[k,k])
```

where z_{α/2} is the standard normal quantile.

**Important**: Due to rotational non-identifiability, we must align the estimated positions X̂ with true positions X_0 via Procrustes rotation before computing coverage.

## Implementation Details

### 1. Data Generation (`generate_latent_curve`)

The true latent positions follow a 3D parametric curve:

```r
t = i/n for i in {1,...,n}
X[i,1] = 0.15 * sin(2πt) + 0.6
X[i,2] = 0.15 * cos(2πt) + 0.6
X[i,3] = 0.15 * cos(4πt)
```

This creates a smooth manifold structure in latent space.

### 2. Covariance Computation (`compute_G_in`)

For vertex i:

1. Compute inner products: s_ij = x_{0i}^T y_{0j} for all j
2. Compute weights: w_j = Ψ'(s_ij) using the dpsi function
3. Network contribution: G_net = Σ_j w_j * (y_{0j} y_{0j}^T)
4. Covariate contribution: G_cov = Z_0^T Z_0
5. Combine: G_in = (G_net + G_cov) / (n + p_cov)

The asymptotic variance is then Σ_i = G_in^{-1}.

### 3. Procrustes Alignment

Since the GRDPG model is only identifiable up to orthogonal rotation, we align estimates to true positions:

1. Compute M = X̂^T X_0
2. SVD: M = U Σ V^T
3. Optimal rotation: R = U V^T
4. Aligned estimates: X_aligned = X̂ R

This alignment is crucial for meaningful coverage testing.

### 4. Coverage Testing

For each replicate:
1. Generate data from the model
2. Fit using fit_grdpg_cov()
3. Align estimates via Procrustes
4. For each test vertex i:
   - Compute G_in using true parameters
   - Invert to get Σ_i
   - Compute standard errors: SE[k] = √(Σ_i[k,k])
   - Check if x_{0i}[k] ∈ [x̂_i[k] ± z_{α/2} * SE[k]]

Empirical coverage is the proportion of replicates where the true value falls in the CI.

## Key Considerations

### 1. Using True Parameters

Notice that G_in is computed using **true** parameters (X_0, Y_0, Z_0), not estimates. This is because:
- We're testing the asymptotic theory, which uses the true parameter values
- In practice, one would plug in estimates, but for simulation validation we use truth

### 2. Scaling Factor

The asymptotic variance includes a √(n+p) scaling. In the covariance computation, we divide by (n+p_cov) to account for this.

### 3. Numerical Stability

We check if G_in is invertible by examining eigenvalues. If min(λ) < 10^{-10}, we skip that vertex to avoid numerical issues.

### 4. Convergence

The Fisher scoring algorithm may not always converge in `maxit` iterations. We track convergence rate across replicates.

## Expected Results

### Coverage Probability

For a nominal 95% confidence level (α = 0.05), we expect:
- Empirical coverage ≈ 95% for each coordinate
- Some variation due to finite samples and Monte Carlo error

If coverage is significantly different from 95%, it suggests:
- Sample size is too small (asymptotic approximation not valid)
- Model misspecification
- Numerical issues in estimation

### Asymptotic Normality

The scatter plots should show:
- Points clustered around the true value (red cross)
- Approximately elliptical shape (reflecting the covariance structure)
- 95% confidence ellipse containing ~95% of points

The Q-Q plots should be approximately linear if the asymptotic normality holds.

## Running the Simulations

### Quick Test
```r
source("simulations/test_coverage_quick.R")
```

This runs a single replicate to verify:
- Data generation works
- Model fitting converges
- G_in is invertible
- Basic coverage check succeeds

### Full Coverage Simulation
```r
source("simulations/coverage_simulation.R")
```

This runs 200 replicates and reports:
- Empirical coverage rates by coordinate
- Overall coverage
- RMSE and standard errors
- Convergence rate

Also tests different sample sizes (n = 40, 60, 80, 100) to see how coverage improves with n.

### Visualization
```r
source("simulations/asymptotic_normality_viz.R")
```

This creates:
- Scatter plots (2D projections of the distribution)
- Marginal histograms with theoretical normal overlay
- Q-Q plots for normality checking

## Interpreting Results

### Good Results
- Coverage rates: 93-97% (allowing for Monte Carlo variation)
- RMSE decreases with increasing n
- Q-Q plots are linear
- Scatter plots show elliptical clusters

### Problematic Results
- Coverage << 95%: Sample size too small or model issues
- Coverage >> 95%: Standard errors may be overestimated
- Non-linear Q-Q plots: Asymptotic approximation not valid
- Irregular scatter plots: Potential issues with alignment or estimation

## Files Generated

After running the simulations:
- `simulations/coverage_results.rds` - Main coverage results
- `simulations/coverage_by_size.rds` - Results for different n
- `simulations/asymptotic_normality_data.rds` - Data for visualization
- `simulations/plots/scatter_asymptotic_normality.pdf` - Scatter plots
- `simulations/plots/histograms_asymptotic_normality.pdf` - Histograms
- `simulations/plots/qq_asymptotic_normality.pdf` - Q-Q plots

## Next Steps

After running the quick test successfully:
1. Run the full coverage simulation
2. Check if empirical coverage ≈ 95%
3. Run the visualization script
4. Examine plots for normality
5. Try different parameter settings if needed

## Potential Extensions

1. **Larger dimensions**: Test with d > 3
2. **Different curves**: Try other parametric curves for X_0
3. **Varying sparsity**: Change the curve amplitude to affect edge probabilities
4. **Covariate strength**: Vary the noise level in B = ZX^T + noise
5. **Plug-in estimates**: Use estimated G_in instead of true values (more realistic)
