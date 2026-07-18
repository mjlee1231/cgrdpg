% Compare three optimization approaches:
% 1. Current: Batch surrogate with trust-region fminunc
% 2. Alternative: Batch surrogate with quasi-newton fminunc
% 3. New: Vectorized Jacobi coordinate ascent

clear; clc;
addpath('core');

fprintf('Comparing Optimization Approaches for GRDPG\n');
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

%% Method 1: Batch Surrogate + Trust-Region
fprintf('======================================\n');
fprintf('Method 1: Batch Surrogate (Trust-Region)\n');
fprintf('======================================\n');

options_tr = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

t0 = tic;
[X_tr, ~, ~, ~, output_tr, ~] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options_tr);
time_tr = toc(t0);

[X_tr_aligned, ~] = procrustes_align(X_tr, X0);
sse_tr = sum((X_tr_aligned - X0).^2, 'all');

fprintf('\nResults: SSE=%.4f, Time=%.1fs, Iters=%d\n\n', ...
    sse_tr, time_tr, output_tr.iterations);

%% Method 2: Batch Surrogate + Quasi-Newton
fprintf('======================================\n');
fprintf('Method 2: Batch Surrogate (Quasi-Newton)\n');
fprintf('======================================\n');

options_qn = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true);

t0 = tic;
[X_qn, ~, ~, ~, output_qn, ~] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options_qn);
time_qn = toc(t0);

[X_qn_aligned, ~] = procrustes_align(X_qn, X0);
sse_qn = sum((X_qn_aligned - X0).^2, 'all');

fprintf('\nResults: SSE=%.4f, Time=%.1fs, Iters=%d\n\n', ...
    sse_qn, time_qn, output_qn.iterations);

%% Method 3: Vectorized Jacobi Coordinate Ascent
fprintf('======================================\n');
fprintf('Method 3: Vectorized Jacobi Coordinate Ascent\n');
fprintf('======================================\n');

t0 = tic;
[X_jacobi, ~, ~, ~, output_jacobi, ~] = ...
    fit_grdpg_jacobi(A, B, d, p, tau, maxit, tol);
time_jacobi = toc(t0);

[X_jacobi_aligned, ~] = procrustes_align(X_jacobi, X0);
sse_jacobi = sum((X_jacobi_aligned - X0).^2, 'all');

fprintf('\nResults: SSE=%.4f, Time=%.1fs, Iters=%d\n\n', ...
    sse_jacobi, time_jacobi, output_jacobi.iterations);

%% Summary Comparison
fprintf('========================================\n');
fprintf('Summary Comparison\n');
fprintf('========================================\n\n');

fprintf('%-30s | %10s | %10s | %10s\n', 'Method', 'SSE', 'Time (s)', 'Iters');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-30s | %10.4f | %10.1f | %10d\n', ...
    'Batch + Trust-Region', sse_tr, time_tr, output_tr.iterations);
fprintf('%-30s | %10.4f | %10.1f | %10d\n', ...
    'Batch + Quasi-Newton', sse_qn, time_qn, output_qn.iterations);
fprintf('%-30s | %10.4f | %10.1f | %10d\n', ...
    'Jacobi Coordinate Ascent', sse_jacobi, time_jacobi, output_jacobi.iterations);
fprintf('%s\n', repmat('-', 1, 70));

% Find best SSE
[best_sse, best_idx] = min([sse_tr, sse_qn, sse_jacobi]);
method_names = {'Batch + Trust-Region', 'Batch + Quasi-Newton', 'Jacobi Coordinate Ascent'};

fprintf('\nBest SSE: %.4f (%s)\n', best_sse, method_names{best_idx});
fprintf('Target (R with Fisher scoring): SSE = 8.23\n');

% Improvements
if sse_jacobi < sse_tr
    fprintf('\nJacobi improvement over Trust-Region: %.1f%%\n', ...
        100 * (sse_tr - sse_jacobi) / sse_tr);
end
