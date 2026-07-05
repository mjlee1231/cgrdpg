% Test CGRDPG with fminunc on 3D simulation
% Compare with R implementation results

clear; clc;

%% Parameters (matching R simulations)
n = 1000;  % number of nodes
p_cov = 1000;  % number of covariates
d = 3;  % embedding dimension
p = 2;  % positive signatures (q = 1 negative)
tau = 0.001;

fprintf('========================================\n');
fprintf('3D CGRDPG Simulation with fminunc\n');
fprintf('========================================\n');
fprintf('n = %d, p_cov = %d, d = %d, tau = %.6f\n\n', n, p_cov, d, tau);

%% Generate synthetic data (same structure as R)
rng(123);  % For reproducibility

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

fprintf('Data generated:\n');
fprintf('  Network density: %.4f\n', sum(A(:)) / (n * (n-1)));
fprintf('  Number of edges: %d\n\n', sum(A(:)) / 2);

%% Fit CGRDPG with fminunc

% Configure fminunc options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true);

% Run optimization
tic;
[X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc(A, B, d, p, tau, options);
elapsed_time = toc;

fprintf('\n========================================\n');
fprintf('Optimization Results\n');
fprintf('========================================\n');
fprintf('Elapsed time: %.2f seconds\n', elapsed_time);
fprintf('Exit flag: %d\n', exitflag);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Final objective: %.6f\n', fval);

%% Compute SSE
sse_matlab = sum((X_opt(:) - X_true(:)).^2);
fprintf('\nSum of Squared Errors (SSE): %.4f\n', sse_matlab);

%% Compare with ASE
fprintf('\n========================================\n');
fprintf('Comparison with ASE\n');
fprintf('========================================\n');

% ASE
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

[V, D] = eig(A_aug);
[eigvals, idx] = sort(diag(D), 'descend');
V = V(:, idx);
X_ase = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

% Align with Procrustes
[~, X_ase_aligned] = procrustes(X_true, X_ase);
[~, X_opt_aligned] = procrustes(X_true, X_opt);

sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);
sse_opt = sum((X_opt_aligned(:) - X_true(:)).^2);

fprintf('SSE (ASE):    %.4f\n', sse_ase);
fprintf('SSE (fminunc): %.4f\n', sse_opt);
fprintf('Improvement:  %.2f%%\n', 100 * (sse_ase - sse_opt) / sse_ase);

%% Save results
results = struct();
results.X_true = X_true;
results.X_opt = X_opt;
results.X_ase = X_ase;
results.Z_true = Z_true;
results.Z_opt = Z_opt;
results.sse_matlab = sse_matlab;
results.sse_ase = sse_ase;
results.sse_opt = sse_opt;
results.fval = fval;
results.exitflag = exitflag;
results.output = output;
results.elapsed_time = elapsed_time;
results.params = struct('n', n, 'p_cov', p_cov, 'd', d, 'p', p, 'tau', tau);

save('matlab_3d_test_results.mat', 'results');
fprintf('\nResults saved to: matlab_3d_test_results.mat\n');
