% Quick test of surrogate likelihood implementation
% Smaller problem for fast local testing

clear; clc;

fprintf('Quick Test: Surrogate Likelihood Implementation\n');
fprintf('================================================\n\n');

%% Small test problem
n = 500;
p_cov = 100;
d = 3;
p = 2;
tau = 0.001;
rng(42);

fprintf('Parameters: n=%d, p_cov=%d, d=%d\n\n', n, p_cov, d);

%% Generate data
X_true = randn(n, d) * 0.5;
Z_true = randn(p_cov, d) * 0.3;
S = diag([1, 1, -1]);

Y_true = X_true * S;
P_net = max(min(Y_true * Y_true', 0.99), 0.01);

A = double(rand(n) < P_net);
A = triu(A, 1);
A = A + A';

B = Z_true * X_true' + randn(p_cov, n) * 0.1;

%% ASE baseline
A_aug = A;
A_aug(1:n+1:end) = sum(A, 2) / (n - 1);  % Set diagonal using linear indexing
[V, D] = eig(A_aug);
[eigvals, idx] = sort(diag(D), 'descend');
X_ase = V(:, idx(1:d)) * diag(sqrt(abs(eigvals(1:d))));

[~, X_ase_aligned] = procrustes(X_true, X_ase);
sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);

fprintf('ASE SSE: %.4f\n\n', sse_ase);

%% Test surrogate
fprintf('Testing surrogate likelihood...\n');
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
    'Display', 'off', 'MaxIterations', 50, ...
    'SpecifyObjectiveGradient', true);

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
