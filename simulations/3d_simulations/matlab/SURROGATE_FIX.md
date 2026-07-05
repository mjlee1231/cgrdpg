# MATLAB Implementation Issue and Fix

## Problem Discovered

The initial MATLAB implementation (`fit_grdpg_fminunc.m`) produced **identical SSE for ASE and fminunc**, converging in just 1 iteration without any improvement.

### Results Comparison:

**R Implementation (n=1000):**
- ASE SSE: ~70-77
- CGRDPG (Fisher) SSE: ~6 ← **90% improvement!**
- OSE SSE: ~200-700 (unstable)

**MATLAB Initial Implementation (n=1000):**
- ASE SSE: ~251
- fminunc SSE: ~251 ← **No improvement!**
- Converged in 1 iteration

## Root Cause

The R `cgrdpg` package uses a **surrogate likelihood** approach, NOT simple probability clamping!

### R's Surrogate Objective:
```r
net <- sum((A - S) * psi(S, tau) + Psi(S, tau))
cov <- -0.5 * sum((B - Z %*% t(X))^2)
objective = net + cov
```

Where:
- `psi(x) = psi1(x) - psi2(x)` is a smoothed derivative
- `Psi(x) = Psi1(x) - Psi2(x)` is the surrogate log-likelihood
- These functions extend the Bernoulli log-likelihood **outside [0,1]** using quadratic patches
- This provides a "barrier-less, twice-differentiable surrogate" that doesn't require constraints

### MATLAB's Original (Incorrect) Objective:
```matlab
nll_network = -sum(A_vec .* log(P_vec) + (1 - A_vec) .* log(1 - P_vec));
sse_cov = sum((B(:) - B_pred(:)).^2);
f = nll_network + 0.5 * sse_cov;
```

This uses:
- Simple **clamping** of probabilities to [tau, 1-tau]
- Standard negative log-likelihood
- Different weighting of covariate term

## Solution

Created new implementation: `fit_grdpg_fminunc_surrogate.m`

**Key components:**

1. **psi_functions.m**: MATLAB versions of R's psi/Psi functions
   - `psi(x)`: Smoothed log(x/(1-x)) with quadratic extensions
   - `Psi(x)`: Integral of psi (the surrogate log-likelihood itself)
   - `dpsi(x)`: Derivative for gradient computation

2. **fit_grdpg_fminunc_surrogate.m**: Corrected optimization
   - Uses same surrogate objective as R
   - Properly handles gradients with respect to the surrogate

3. **test_surrogate_vs_clamped.m**: Comparison script
   - Tests both implementations side-by-side
   - Should show ~90% improvement with surrogate approach

## Testing

Before re-running 100 replications:

```matlab
% Run local test
cd simulations/3d_simulations/matlab
matlab -nodisplay -nosplash -r "test_surrogate_vs_clamped; exit;"
```

Expected outcome:
- Surrogate implementation should achieve SSE ~6-10 (similar to R)
- Clamped implementation will stay at SSE ~250 (no improvement)

## Next Steps

1. ✅ Test `test_surrogate_vs_clamped.m` locally
2. ⬜ Update 100-replication scripts to use `fit_grdpg_fminunc_surrogate`
3. ⬜ Re-run HPC jobs with corrected implementation
4. ⬜ Compare MATLAB surrogate results with R Fisher results

## References

- R implementation: `R/psi.R` - psi/Psi functions
- R implementation: `R/fisher.R` - surrogate objective and Fisher scoring
- R implementation: `R/fit.R` - main fitting function using surrogate
