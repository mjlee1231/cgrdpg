# Parallel Optimization Summary

## Overview
Optimized model fitting (`fit_grdpg_cov`) using parallel processing with `foreach` and `doParallel` packages.

## Benchmark Results (n=200, 10 iterations)

| Metric | Serial | Parallel (11 cores) | Improvement |
|--------|--------|---------------------|-------------|
| **Time** | 44.18 sec | 11.53 sec | **3.83x faster** |
| **Parallel Efficiency** | - | 34.8% | - |
| **SSE** | 71.39 | 70.25 | Similar convergence |

## Key Findings

### 1. Significant Speedup
- **3.83x faster** with 11 cores on Apple M2 Max
- Parallel efficiency: 34.8% (reasonable for this workload)
- Speedup will scale with more cores available on HPC

### 2. Extrapolation to Full Simulation (n=1000, 100 reps)

**Per replication:**
- Serial: ~55.2 minutes
- Parallel: ~14.4 minutes
- **Time saved: ~40.8 minutes per replication**

**For 100 replications:**
- **Total time saved: ~68 hours (2.8 days)**
- From ~92 hours → ~24 hours with parallelization

### 3. Implementation Details

**Parallelization Strategy:**
- Parallel sweep uses **Jacobi-style updates** (all nodes updated simultaneously)
- Serial version uses **Gauss-Seidel** (sequential updates with immediate incorporation)
- Difference: Jacobi may require slightly more iterations but is much faster overall

**Bottleneck Parallelized:**
- `fisher_sweep_X` function's `for (i in 1:n)` loop
- Each node's Fisher scoring update computed independently
- Uses `foreach` with `.combine = rbind` to aggregate results

## Files Created

1. **R/fisher_parallel.R** - Parallel implementations
   - `fisher_sweep_X_parallel()` - Parallelized Fisher sweep
   - `fit_grdpg_cov_parallel()` - Parallelized model fitting

2. **simulations/benchmark_parallel_fitting.R** - Comprehensive benchmark script

## Usage

```r
library(cgrdpg)

# Parallel version (recommended for large n)
fit <- fit_grdpg_cov_parallel(
  A, B, d = 2, p = 2, q = 0,
  tau = 0.005, maxit = 30, tol = 0.01,
  ncores = 11  # or NULL to auto-detect
)

# Serial version (original)
fit <- fit_grdpg_cov(
  A, B, d = 2, p = 2, q = 0,
  tau = 0.005, maxit = 30, tol = 0.01
)
```

## Recommendations

### For HPC Simulations (n=1000, 100 reps):
✅ **HIGHLY RECOMMENDED** to use parallel version
- Allocate 16+ CPUs per job
- Expected speedup: 4-8x depending on available cores
- Will reduce total computation time from ~4 days to <1 day

### Dependencies:
- `foreach` package
- `doParallel` package
- Both installed automatically on first run

## Next Steps

### Option 1: Use Parallel Version on HPC ✅ RECOMMENDED
- Modify SLURM script to allocate more CPUs (e.g., 16-32)
- Update simulation script to use `fit_grdpg_cov_parallel()`
- Expected total time: ~24 hours for 100 reps (vs 92 hours serial)

### Option 2: Rcpp Optimization (if more speed needed)
- Convert core Fisher scoring loop to C++
- Potential for additional 2-5x speedup
- More complex to implement and maintain

## Technical Notes

**Parallel Efficiency:**
- 34.8% efficiency is reasonable for this workload
- Limited by:
  - Communication overhead between processes
  - Memory bandwidth (large matrices)
  - Synchronization at each iteration

**Convergence:**
- Jacobi parallelization may differ slightly from Gauss-Seidel
- Both converge to similar solutions (SSE within 1.5%)
- Parallel may require 1-2 more iterations but still much faster overall

**Memory:**
- Each worker needs a copy of matrices (A, X, Y, Z)
- For n=1000: ~8MB per worker (manageable)
- Total memory scales linearly with number of cores

## Profiling Hierarchy

From slowest to fastest optimization approach:

1. ❌ **No optimization** - 55 min/rep @ n=1000
2. ✅ **Vectorize G_in computation** - Saves ~30 sec/rep (0.5% improvement)
3. ✅ **Parallel Fisher scoring** - Saves ~41 min/rep (3.83x speedup) ← **CURRENT**
4. 🔧 **Rcpp conversion** - Potential additional 2-5x (future work)

## Conclusion

**Parallelization provides the single biggest performance improvement** with minimal code changes. The 3.83x speedup translates to saving 2.8 days of computation time for the full simulation. Strongly recommended for HPC deployment.
