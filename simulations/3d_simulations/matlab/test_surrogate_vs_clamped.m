% Compare clamped likelihood vs surrogate likelihood implementations
% This test checks if the surrogate approach (matching R) improves over ASE

clear; clc;

fprintf('========================================\n');
fprintf('Testing Surrogate vs Clamped Objective\n');
fprintf('========================================\n\n');

%% Parameters
n = 1000;
p_cov = 1000;
d = 3;
p = 2;  % positive signatures (q = 1 negative)
tau = 0.001;
rng(123);  % For reproducibility

fprintf('Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n\n', n, p_cov, d, tau);

%% Generate synthetic data
fprintf('Generating synthetic data...\n');

% True latent positions (3D: 2 positive, 1 negative)
X_true = randn(n, d) * 0.5;

% True covariate coefficients
Z_true = randn(p_cov, d) * 0.3;

% Signature matrix
S = diag([1, 1, -1]);

% Generate probabilities
Y_true = X_true * S;
P_net = Y_true * Y_true';
P_net = max(min(P_net, 0.99), 0.01);  % Clamp to valid range

% Generate adjacency matrix
A = double(rand(n) < P_net);
A = triu(A, 1);  % Upper triangle only (undirected)
A = A + A';  % Make symmetric

% Generate covariates
B = Z_true * X_true' + randn(p_cov, n) * 0.1;

fprintf('Data generated successfully.\n\n');

%% Compute ASE baseline
fprintf('Computing ASE baseline...\n');
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

[V, D] = eig(A_aug);
[eigvals, idx] = sort(diag(D), 'descend');
V = V(:, idx);
X_ase = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

% Procrustes alignment for ASE
[~, X_ase_aligned] = procrustes(X_true, X_ase);
sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);

fprintf('ASE SSE: %.4f\n\n', sse_ase);

%% Configure fminunc options (suppress iteration display for cleaner output)
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'off', ...  % Suppress iteration output
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true);

%% Test 1: Clamped likelihood (original implementation)
fprintf('========================================\n');
fprintf('Test 1: Clamped Likelihood\n');
fprintf('========================================\n');
t0 = tic;
[X_clamped, Z_clamped, fval_clamped, exitflag_clamped, output_clamped] = ...
    fit_grdpg_fminunc(A, B, d, p, tau, options);
time_clamped = toc(t0);

[~, X_clamped_aligned] = procrustes(X_true, X_clamped);
sse_clamped = sum((X_clamped_aligned(:) - X_true(:)).^2);

fprintf('Results:\n');
fprintf('  SSE: %.4f\n', sse_clamped);
fprintf('  Improvement over ASE: %.2f%%\n', 100 * (sse_ase - sse_clamped) / sse_ase);
fprintf('  Converged: %d\n', exitflag_clamped > 0);
fprintf('  Iterations: %d\n', output_clamped.iterations);
fprintf('  Time: %.2f seconds\n\n', time_clamped);

%% Test 2: Surrogate likelihood (matching R)
fprintf('========================================\n');
fprintf('Test 2: Surrogate Likelihood (R-style)\n');
fprintf('========================================\n');
t0 = tic;
[X_surrogate, Z_surrogate, fval_surrogate, exitflag_surrogate, output_surrogate] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);
time_surrogate = toc(t0);

[~, X_surrogate_aligned] = procrustes(X_true, X_surrogate);
sse_surrogate = sum((X_surrogate_aligned(:) - X_true(:)).^2);

fprintf('Results:\n');
fprintf('  SSE: %.4f\n', sse_surrogate);
fprintf('  Improvement over ASE: %.2f%%\n', 100 * (sse_ase - sse_surrogate) / sse_ase);
fprintf('  Converged: %d\n', exitflag_surrogate > 0);
fprintf('  Iterations: %d\n', output_surrogate.iterations);
fprintf('  Time: %.2f seconds\n\n', time_surrogate);

%% Summary comparison
fprintf('========================================\n');
fprintf('Summary Comparison\n');
fprintf('========================================\n');
fprintf('Method          | SSE    | Improvement | Iterations\n');
fprintf('----------------|--------|-------------|------------\n');
fprintf('ASE             | %.4f | --          | --\n', sse_ase);
fprintf('Clamped         | %.4f | %.1f%%       | %d\n', ...
    sse_clamped, 100 * (sse_ase - sse_clamped) / sse_ase, output_clamped.iterations);
fprintf('Surrogate (R)   | %.4f | %.1f%%       | %d\n', ...
    sse_surrogate, 100 * (sse_ase - sse_surrogate) / sse_ase, output_surrogate.iterations);
fprintf('\n');

if sse_surrogate < sse_ase
    fprintf('SUCCESS: Surrogate approach improves over ASE!\n');
    fprintf('This matches the R implementation behavior.\n');
else
    fprintf('WARNING: Surrogate approach did not improve over ASE.\n');
    fprintf('Check gradient computation or optimization settings.\n');
end

fprintf('\n========================================\n');
fprintf('Expected R-style improvement: ~90%% (SSE ~7 from ~70)\n');
fprintf('========================================\n');
