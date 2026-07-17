% Test Fisher scoring implementation on one replication
% Compare SSE and coverage with expected R results

clear; clc;
addpath('core');

fprintf('Testing Fisher Scoring Implementation\n');
fprintf('======================================\n\n');

% Parameters (match R)
n = 1000;
p_cov = 500;
d = 3;
p = 2;
q = 1;
maxit = 30;
tol = 0.005;  % Match R default
tau = 0.001;
eps_clip = 1e-10;
chi2_crit = chi2inv(0.95, d);

% Set seed (rep 1)
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
fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

% Generate network and covariates
A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';
A(1:n+1:end) = 0;

B = Z0 * X0' + randn(p_cov, n);

% Fit with Fisher scoring
fprintf('Fitting with Fisher scoring...\n');
fprintf('----------------------------------------\n');
t_start = tic;
[X_opt, Z_opt, fval, exitflag, output, S_est] = ...
    fit_grdpg_fisher(A, B, d, p, tau, maxit, tol);
fit_time = toc(t_start);
fprintf('----------------------------------------\n');

% Procrustes alignment
[X_fisher, ~] = procrustes_align(X_opt, X0);
sse_fisher = sum((X_fisher - X0).^2, 'all');

fprintf('\n');
fprintf('Results:\n');
fprintf('  SSE: %.4f (R rep 1: 8.23)\n', sse_fisher);
fprintf('  Time: %.1f sec\n', fit_time);
fprintf('  Converged: %d\n', output.converged);
fprintf('  Iterations: %d\n', output.iterations);

% Compute coverage for a few vertices
fprintf('\nComputing coverage for 10 sample vertices...\n');
Y_fisher = X_fisher * S_est;
Z_fisher = (X_fisher \ B')';

sample_vertices = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900];
n_covered = 0;

for idx = 1:length(sample_vertices)
    i = sample_vertices(idx);
    err = X0(i,:) - X_fisher(i,:);
    G_true = compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau);
    covered = check_coverage(err, G_true, chi2_crit, n + p_cov);
    if covered
        n_covered = n_covered + 1;
    end
end

fprintf('  Sample coverage: %d/10 = %.0f%%\n', n_covered, 100 * n_covered / 10);
fprintf('\nExpected from R (rep 1):\n');
fprintf('  SSE: 8.23\n');
fprintf('  Coverage: ~77%%\n');

if sse_fisher < 15
    fprintf('\n✓ SUCCESS: SSE looks good!\n');
else
    fprintf('\n✗ WARNING: SSE still high, may need debugging\n');
end
