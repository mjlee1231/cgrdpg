% Compare MATLAB fminunc results with R Fisher scoring results
% Reads .rds files from R simulations and compares SSE distributions

clear; clc;

fprintf('========================================\n');
fprintf('Comparing MATLAB fminunc vs R Fisher\n');
fprintf('========================================\n\n');

%% Read R results
% You need to export R results to .mat format first
% In R, run:
%   library(R.matlab)
%   results <- readRDS("results_3d_ase_ose_cgrdpg_n1000/rep_001.rds")
%   writeMat("r_results_rep001.mat", results=results)

r_results_file = '../simulations/3d_simulations/r_results_rep001.mat';

if exist(r_results_file, 'file')
    fprintf('Loading R results from: %s\n', r_results_file);
    r_data = load(r_results_file);

    % Extract SSE values
    sse_r_ase = r_data.results.sse.ase;
    sse_r_ose = r_data.results.sse.ose;
    sse_r_cgrdpg = r_data.results.sse.cgrdpg;

    fprintf('R Results (single replication):\n');
    fprintf('  SSE (ASE):    %.4f\n', sse_r_ase);
    fprintf('  SSE (OSE):    %.4f\n', sse_r_ose);
    fprintf('  SSE (CGRDPG): %.4f\n', sse_r_cgrdpg);
else
    fprintf('R results file not found. Skipping comparison.\n');
    fprintf('To create R results .mat file:\n');
    fprintf('  1. In R: library(R.matlab)\n');
    fprintf('  2. results <- readRDS("path/to/rep_001.rds")\n');
    fprintf('  3. writeMat("r_results_rep001.mat", results=results)\n');
end

%% Load MATLAB results
if exist('matlab_3d_test_results.mat', 'file')
    fprintf('\nLoading MATLAB results...\n');
    matlab_data = load('matlab_3d_test_results.mat');

    fprintf('MATLAB Results (fminunc):\n');
    fprintf('  SSE (ASE):     %.4f\n', matlab_data.results.sse_ase);
    fprintf('  SSE (fminunc): %.4f\n', matlab_data.results.sse_opt);
    fprintf('  Optimization time: %.2f seconds\n', matlab_data.results.elapsed_time);
    fprintf('  Exit flag: %d\n', matlab_data.results.exitflag);
    fprintf('  Iterations: %d\n', matlab_data.results.output.iterations);
else
    fprintf('\nMATLAB results not found. Run test_3d_simulation.m first.\n');
end

%% Summary
fprintf('\n========================================\n');
fprintf('Summary\n');
fprintf('========================================\n\n');

fprintf('Key advantages of MATLAB fminunc:\n');
fprintf('  1. Multiple optimization algorithms (quasi-Newton, trust-region)\n');
fprintf('  2. Automatic gradient checking\n');
fprintf('  3. Robust line search with strong Wolfe conditions\n');
fprintf('  4. Better numerical conditioning\n');
fprintf('  5. Detailed convergence diagnostics\n\n');

fprintf('Next steps:\n');
fprintf('  1. Run multiple replications with different random seeds\n');
fprintf('  2. Compare SSE distributions (MATLAB vs R)\n');
fprintf('  3. Check for numerical instability (NaN, Inf values)\n');
fprintf('  4. Try different fminunc algorithms:\n');
fprintf('     - quasi-newton (default, good for large problems)\n');
fprintf('     - trust-region (requires Hessian, more stable)\n');
