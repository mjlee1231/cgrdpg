% Quick test of surrogate likelihood implementation
% Smaller problem for fast local testing

clear; clc;

% Add core folder to path
addpath('core');

fprintf('Quick Test: Surrogate Likelihood Implementation\n');
fprintf('================================================\n\n');

%% Parameters matching R simulation (ase_ose_cgrdpg_vertex_wise_3d_n1000.R)
n = 500;
p_cov = 250;  % Using 250 for quick test (R uses 500 for n=1000)
d = 3;
p = 2;
tau = 0.001;
rng(598);  % Same seed as R

fprintf('Parameters: n=%d, p_cov=%d, d=%d\n', n, p_cov, d);
fprintf('(Matching R simulation data generation)\n\n');

%% Generate data EXACTLY like R simulation
% Latent positions: specific pattern (not random!)
t = (1:n)' / n;
X_true = [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)];

% Signature matrix
S = diag([1, 1, -1]);

% Covariate coefficients
Z_true = randn(p_cov, d);

% Edge probabilities
Y_true = X_true * S;
P_net = X_true * Y_true';

fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P_net(:)), max(P_net(:)));

% Generate adjacency matrix (undirected)
A = double(rand(n) < P_net);
A = triu(A, 1);  % Upper triangle only
A = A + A';  % Make symmetric

% Generate covariates with noise sd=1.0 (matching R)
B = Z_true * X_true' + randn(p_cov, n);

%% ASE baseline
A_aug = A;
A_aug(1:n+1:end) = sum(A, 2) / (n - 1);  % Set diagonal using linear indexing
[V, D] = eig(A_aug);
eigvals = diag(D);

% Sort by MAGNITUDE (matching R's eigs_sym with which="LM")
[~, idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(idx);
V = V(:, idx);

X_ase = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

[~, X_ase_aligned] = procrustes(X_true, X_ase);
sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);

fprintf('ASE SSE: %.4f\n\n', sse_ase);

%% Test surrogate
fprintf('Testing surrogate likelihood...\n');

% R uses maxit=30, tol=0.01 (max row change)
% Use trust-region algorithm for better numerical stability
options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 30, ...  % Match R's maxit
    'OptimalityTolerance', 1e-4, ...  % Loosen from default 1e-6
    'StepTolerance', 1e-8, ...  % Loosen from default 1e-10
    'FunctionTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'lbfgs');

tic;
[X_opt, ~, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);
t_elapsed = toc;

[~, X_opt_aligned] = procrustes(X_true, X_opt);
sse_opt = sum((X_opt_aligned(:) - X_true(:)).^2);

fprintf('\nResults:\n');
fprintf('  Surrogate SSE: %.4f\n', sse_opt);
fprintf('  Improvement: %.1f%%\n', 100 * (sse_ase - sse_opt) / sse_ase);
fprintf('  Iterations: %d\n', output.iterations);
fprintf('  Time: %.2f sec\n', t_elapsed);
fprintf('  Exit flag: %d\n\n', exitflag);

if sse_opt < 0.5 * sse_ase
    fprintf('✓ SUCCESS: Surrogate improved over ASE by >50%%\n');
else
    fprintf('✗ WARNING: Limited improvement (expected ~90%%)\n');
end
