% Test ONLY Jacobi coordinate ascent method
% Quick test to verify it works without numerical explosion

clear; clc;
addpath('core');

fprintf('Testing Vectorized Jacobi Coordinate Ascent\n');
fprintf('===========================================\n\n');

% Parameters (rep 1)
n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;
maxit = 30;
tol = 0.005;

% Set seed
rng(598 + 1);

% Generate data
t = (1:n)' / n;
X0 = [0.15 * sin(2*pi*t) + 0.6, ...
      0.15 * cos(2*pi*t) + 0.6, ...
      0.15 * cos(4*pi*t)];

S = diag([1, 1, -1]);
Y0 = X0 * S;
Z0 = randn(p_cov, d);
P = X0 * Y0';

A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';
A(1:n+1:end) = 0;
B = Z0 * X0' + randn(p_cov, n);

fprintf('Data: n=%d, p_cov=%d, d=%d\n', n, p_cov, d);
fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

%% Jacobi Method
t0 = tic;
[X_jacobi, ~, ~, ~, output_jacobi, ~] = ...
    fit_grdpg_jacobi(A, B, d, p, tau, maxit, tol);
time_jacobi = toc(t0);

[X_jacobi_aligned, ~] = procrustes_align(X_jacobi, X0);
sse_jacobi = sum((X_jacobi_aligned - X0).^2, 'all');

fprintf('\n========================================\n');
fprintf('Results\n');
fprintf('========================================\n');
fprintf('SSE:        %.4f\n', sse_jacobi);
fprintf('Time:       %.1f sec\n', time_jacobi);
fprintf('Iterations: %d\n', output_jacobi.iterations);
fprintf('Converged:  %s\n', mat2str(output_jacobi.converged));

fprintf('\n========================================\n');
fprintf('Comparison with Baselines\n');
fprintf('========================================\n');
fprintf('Target (R Fisher scoring):    SSE = 8.23\n');
fprintf('Batch Trust-Region (MATLAB):  SSE = 51.15\n');
fprintf('Jacobi (this run):            SSE = %.2f\n', sse_jacobi);

if sse_jacobi < 51.15
    improvement = 100 * (51.15 - sse_jacobi) / 51.15;
    fprintf('\nJacobi improvement over batch: %.1f%%\n', improvement);
end

if sse_jacobi > 100
    fprintf('\n⚠ WARNING: SSE > 100 suggests numerical issues\n');
elseif sse_jacobi < 30
    fprintf('\n✓ SUCCESS: SSE < 30 is competitive with R!\n');
elseif sse_jacobi < 51.15
    fprintf('\n✓ GOOD: SSE better than batch surrogate\n');
else
    fprintf('\n✗ Poor: SSE worse than batch surrogate\n');
end
