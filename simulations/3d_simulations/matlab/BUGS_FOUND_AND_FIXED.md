# Critical Bugs Found and Fixed in MATLAB Implementation

## Summary

The initial MATLAB implementation had **3 critical bugs** that prevented it from matching R's performance. After fixes, we expect ~90% improvement over ASE (matching R).

---

## Bug 1: Wrong Objective Function (Most Critical!)

**Impact:** Complete failure - 0% improvement over ASE

**Problem:**
Original `fit_grdpg_fminunc.m` used **probability clamping** instead of R's **surrogate likelihood**.

**R's approach:**
```r
# Surrogate objective with smooth psi/Psi functions
net <- sum((A - S) * psi(S, tau) + Psi(S, tau))
cov <- -0.5 * sum((B - Z %*% t(X))^2)
objective <- net + cov
```

**MATLAB's wrong approach:**
```matlab
% Simple clamping - NOT a surrogate!
P_clamped = max(min(P_net, 1 - tau), tau);
nll = -sum(A .* log(P_clamped) + (1-A) .* log(1-P_clamped));
```

**Fix:**
Created `fit_grdpg_fminunc_surrogate.m` with proper `psi_functions.m` that implements:
- `psi(x)`: Smoothed log(x/(1-x)) with quadratic extensions outside [0,1]
- `Psi(x)`: Integral of psi (the surrogate log-likelihood itself)
- `dpsi(x)`: Derivative for gradient computation

---

## Bug 2: Eigenvalue Sorting by Value Instead of Magnitude

**Impact:** Wrong ASE initialization, especially harmful for GRDPG with negative eigenvalues

**Problem:**
MATLAB sorted eigenvalues by **descending value**, but R sorts by **descending magnitude**.

**Example:**
Eigenvalues: [5, -4, 3, -2, 1]
- **R (magnitude):** Selects [5, -4, 3] ✓ Keeps important negative eigenvalue
- **MATLAB (value):** Selects [5, 3, 1] ✗ Misses the -4 eigenvalue!

**R's approach:**
```r
eigenA <- RSpectra::eigs_sym(A, k = d, which = "LM")  # "LM" = largest magnitude
```

**MATLAB's wrong approach:**
```matlab
[eigvals, idx] = sort(diag(D), 'descend');  # Sorts by value, not magnitude!
```

**Fix:**
```matlab
% Sort by MAGNITUDE (absolute value)
[~, idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(idx);
V = V(:, idx);
```

---

## Bug 3: Objective Summed Over Upper Triangle Instead of Full Matrix

**Impact:** Wrong objective value and gradient (off by factor of 2)

**Problem:**
MATLAB summed network term over **upper triangle only**, but R sums over **all (i,j) pairs with i≠j**.

For undirected graphs, this means:
- **R:** Counts each edge twice (as i→j and j→i)
- **MATLAB:** Counts each edge once

**R's approach:**
```r
S <- X %*% t(Y)
diag(S) <- 0  # Exclude diagonal
net <- sum((A - S) * psi(S) + Psi(S))  # Sums over FULL n×n matrix
```

**MATLAB's wrong approach:**
```matlab
triu_idx = triu(true(n), 1);  # Upper triangle only
net_obj = sum((A_vec - S_vec) .* psi_vec + Psi_vec);  # Half the value!
```

**Fix:**
```matlab
% Sum over full matrix (excluding diagonal) to match R
S_mat(1:n+1:end) = 0;  % Zero diagonal
net_obj = sum((A(:) - S_mat(:)) .* psi_val(:) + Psi_val(:));  # Full matrix
```

**Gradient fix:**
Added factor of 2 because derivative includes both i→j and j→i terms:
```matlab
W_net = (A - S_mat) .* dpsi_val;
W_net(1:n+1:end) = 0;
grad_X_net = -2 * W_net * Y * S;  # Factor of 2 from symmetric summation
```

---

## Files Updated

### Corrected Implementation:
- ✅ `psi_functions.m` - NEW: R's psi/Psi functions
- ✅ `fit_grdpg_fminunc_surrogate.m` - NEW: Corrected surrogate implementation
- ✅ `quick_test_surrogate.m` - NEW: Fast test script

### Need to Update (still using wrong version):
- ❌ `run_100_replications_parallel.m` - n=1000 parallel
- ❌ `run_100_replications_n2000_parallel.m` - n=2000 parallel
- ❌ `run_100_replications.m` - serial version
- ❌ `test_3d_simulation.m` - single test

These files still call `fit_grdpg_fminunc` (wrong) instead of `fit_grdpg_fminunc_surrogate` (correct).

---

## Next Steps

1. ✅ Test `quick_test_surrogate.m` on HPC (expect ~90% improvement)
2. ⬜ Update 100-rep scripts to use `fit_grdpg_fminunc_surrogate`
3. ⬜ Re-run both n=1000 and n=2000 jobs
4. ⬜ Compare MATLAB vs R results (should now match!)

---

## Expected Results

After all fixes:

**R Results (n=1000):**
- ASE SSE: ~70-77
- Fisher SSE: ~6 (90% improvement)

**MATLAB Results (n=1000) - After Fixes:**
- ASE SSE: ~70-77 (should match R now!)
- Surrogate SSE: ~6 (should match R!)
- Improvement: ~90%

**Previous MATLAB (With Bugs):**
- ASE SSE: ~251 (wrong due to eigenvalue sorting bug)
- fminunc SSE: ~251 (no improvement due to wrong objective)
