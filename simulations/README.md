# Latent Position Coverage Probability Simulations

This directory contains simulation scripts for testing the empirical coverage probability of vertex-wise confidence sets for latent position vectors.

## Files

1. **test_coverage_quick.R** - Quick test to verify the setup is viable
2. **coverage_simulation.R** - Full simulation for empirical coverage rates
3. **asymptotic_normality_viz.R** - Visualization of asymptotic normality

## Running the Simulations

### Prerequisites

Install the package first (from the repository root):

```r
# Install dependencies if needed
install.packages(c("MASS", "ggplot2", "reshape2"))

# Install the cgrdpg package
devtools::install()
# Or if in RStudio:
# devtools::load_all()
```

### Quick Test

Run the quick test to verify the setup:

```r
source("simulations/test_coverage_quick.R")
```

This will:
- Generate a single dataset with the latent curve structure
- Fit the GRDPG model
- Compute the covariance matrix G_in^{-1}
- Check confidence interval coverage for one vertex

### Full Coverage Simulation

Run the full coverage probability simulation:

```r
source("simulations/coverage_simulation.R")
```

This performs Monte Carlo simulations to:
- Estimate empirical coverage rates across multiple replicates
- Test different sample sizes and parameter settings
- Generate coverage probability plots

### Asymptotic Normality Visualization

Visualize the asymptotic normality:

```r
source("simulations/asymptotic_normality_viz.R")
```

This creates scatter plots showing the distribution of estimated positions across repeated experiments.

## Key Formulas

The covariance matrix for vertex i is:

```
Σ_i = G_in^{-1}
```

where

```
G_in = (1/(n+p)) * [Σ_{j=1}^n Ψ'(x_{0i}^T y_{0j}) y_{0j} y_{0j}^T + Σ_{l=1}^p z_{0l} z_{0l}^T]
```

The (1-α)% confidence region for x_i (up to rotation) is based on asymptotic normality:

```
√(n+p) * (x̂_i - x_{0i}) → N(0, Σ_i)
```
