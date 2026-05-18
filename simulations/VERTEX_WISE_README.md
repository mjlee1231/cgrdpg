# Vertex-wise Coverage Simulation (n=500, 100 Replications)

This directory contains scripts for computing vertex-wise coverage rates across 100 replications.

## Overview

- **Goal**: Compute coverage rate for each of 500 vertices across 100 independent replications
- **Method**: Elliptical confidence intervals using Oracle Procrustes alignment
- **Configuration**: 2D SPREAD OPTIMIZED (r=0.28, c=0.42)
  - Edge probability range: [0.195, 0.764]
- **Parameters**: n=500, p_cov=250, d=2, tau=0.005, seed=598

## File Structure

```
simulations/
├── vertex_wise_coverage_spread_100reps.R     # Main simulation script (run per replication)
├── submit_vertex_wise_spread_100reps.slurm   # SLURM array job submission script
├── aggregate_vertex_wise_coverage.R          # Aggregates results from all reps
├── vertex_wise_results_n500_optimized/       # Output directory (created automatically)
│   ├── vertex_wise_spread_rep001.rds         # Results from rep 1
│   ├── vertex_wise_spread_rep002.rds         # Results from rep 2
│   ├── ...
│   ├── vertex_wise_spread_rep100.rds         # Results from rep 100
│   ├── vertex_wise_coverage_aggregated.rds   # Final aggregated results
│   └── vertex_wise_coverage_results.csv      # CSV for easy viewing
└── vertex_wise_logs_n500_optimized/          # Log directory (must create manually)
    ├── vertex_wise_spread_rep001_*.out       # Standard output
    ├── vertex_wise_spread_rep001_*.err       # Error output
    └── ...
```

## Usage Instructions

### Step 1: Create Log Directory (On HPC)

```bash
cd ~/cgrdpg/simulations
mkdir -p vertex_wise_logs_n500_optimized
```

### Step 2: Submit Array Job

```bash
sbatch submit_vertex_wise_spread_100reps.slurm
```

This submits 100 parallel jobs (array 1-100), each computing coverage for all 500 vertices.

### Step 3: Monitor Jobs

```bash
# Check job status
squeue -u $USER

# Count running jobs
watch -n 30 'squeue -u $USER | wc -l'

# Check specific job output
tail -f vertex_wise_logs_n500_optimized/vertex_wise_spread_rep001_*.out
```

### Step 4: Aggregate Results

After all 100 jobs complete:

```bash
Rscript aggregate_vertex_wise_coverage.R
```

This creates:
- `vertex_wise_results_n500_optimized/vertex_wise_coverage_aggregated.rds` (full results)
- `vertex_wise_results_n500_optimized/vertex_wise_coverage_results.csv` (summary table)

### Step 5: Transfer Results to Local

```bash
# On local machine
cd ~/Documents/GitHub/cgrdpg/simulations
scp -r mle6@bigred.uits.iu.edu:~/cgrdpg/simulations/vertex_wise_results_n500_optimized .
```

## Expected Runtime

- **Per replication**: ~15-20 minutes (fit + 500 vertex coverage computations)
- **Wall time**: ~20 minutes (all jobs run in parallel on separate nodes)
- **Total CPU time**: ~100 × 20 min = 2000 minutes = 33 hours

## Output Format

### Per-replication file (vertex_wise_spread_repXXX.rds)
```r
list(
  rep = 1,                          # Replication number
  n = 500,
  p_cov = 250,
  d = 2,
  in_ellipse_true = c(T,F,T,...),   # Length 500: coverage indicator for each vertex
  in_ellipse_plugin = c(T,T,F,...), # Length 500: coverage indicator for each vertex
  SSE = 12.34,
  fit_time = 123.4,
  coverage_true = 0.934,            # Overall coverage for this rep
  coverage_plugin = 0.921,
  ...
)
```

### Aggregated file (vertex_wise_coverage_aggregated.rds)
```r
list(
  n = 500,
  n_reps = 100,
  vertex_coverage_true = c(0.95, 0.93, ...),    # Length 500: coverage rate per vertex
  vertex_coverage_plugin = c(0.92, 0.91, ...),  # Length 500: coverage rate per vertex
  in_ellipse_true_matrix = matrix(...),         # 500 × 100 matrix of indicators
  in_ellipse_plugin_matrix = matrix(...),       # 500 × 100 matrix of indicators
  ...
)
```

## Interpretation

**Vertex-wise coverage rate**: For each vertex i, the fraction of replications where vertex i was covered:
```
coverage[i] = (# times vertex i was covered) / 100
```

This allows you to see if certain vertices (e.g., endpoints vs. middle of arc) have systematically different coverage rates.

## Troubleshooting

### Jobs fail with out-of-memory
Increase memory in SLURM script: `#SBATCH --mem=64G`

### Missing result files
Check error logs: `cat vertex_wise_logs_n500/vertex_wise_spread_rep*_*.err`

### Aggregation fails
Make sure all 100 result files exist in `vertex_wise_results_n500/`
