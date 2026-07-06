% Test surrogate likelihood with EXACT specification
% Parameters: n=1000, p_cov=500, d=3, max_iter=30
% Matches R simulation setup exactly

clear; clc;

fprintf('========================================\n');
fprintf('Surrogate Likelihood Test\n');
fprintf('========================================\n\n');

%% Parameters (EXACT specification)
n = 1000;
p_cov = 500;
d = 3;
p = 2;  % Number of positive eigenvalues (q = d - p = 1)
tau = 0.001;
max_iter = 30;
rng(598);  % Fixed seed for reproducibility

fprintf('Parameters:\n');
fprintf('  n = %d (number of nodes)\n', n);
fprintf('  p_cov = %d (number of covariates)\n', p_cov);
fprintf('  d = %d (embedding dimension)\n', d);
fprintf('  p = %d, q = %d (signature: +,+,-)\n', p, d-p);
fprintf('  tau = %.6f (smoothing parameter)\n', tau);
fprintf('  max_iter = %d\n\n', max_iter);

%% Generate data EXACTLY matching R specification
fprintf('Generating data (matching R specification)...\n');

% 1. True latent positions: parametric curve
%    R code: X0 <- cbind(0.15 * sin(2*pi*t) + 0.6,
%                        0.15 * cos(2*pi*t) + 0.6,
%                        0.15 * cos(4*pi*t))
t = (1:n)' / n;
X_true = [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)];

% 2. Signature matrix: S = diag(c(1, 1, -1))
S = diag([1, 1, -1]);

% 3. Covariate coefficients: Z0 <- matrix(rnorm(p_cov * d), p_cov, d)
Z_true = randn(p_cov, d);

% 4. Edge probabilities: Y0 <- X0 %*% S, P <- X0 %*% t(Y0)
Y_true = X_true * S;
P_net = X_true * Y_true';

fprintf('  Edge probability range: [%.4f, %.4f]\n', min(P_net(:)), max(P_net(:)));

% 5. Adjacency matrix: A <- (runif(n^2) < P) * 1.0
%                       A <- A * upper.tri(A, diag = FALSE) + t(...)
A = double(rand(n) < P_net);
A = triu(A, 1);  % Upper triangle only (no diagonal)
A = A + A';      % Make symmetric

fprintf('  Number of edges: %d\n', sum(A(:))/2);

% 6. Covariate matrix: B <- Z0 %*% t(X0) + matrix(rnorm(p_cov * n, sd = 1.0), p_cov, n)
B = Z_true * X_true' + randn(p_cov, n);

fprintf('Data generation complete.\n\n');

%% ASE Baseline (for comparison)
fprintf('Computing ASE baseline...\n');

% Augmented adjacency (diagonal = regularized degree)
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

% Eigendecomposition
[V, D] = eig(A_aug);
eigvals = diag(D);

% Sort by MAGNITUDE (critical for GRDPG!)
[~, idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(idx);
V = V(:, idx);

fprintf('  Top %d eigenvalues (by magnitude):\n', d);
for i = 1:d
    fprintf('    λ_%d = %+.4f (|λ| = %.4f)\n', i, eigvals(i), abs(eigvals(i)));
end

% ASE initialization: X = V_d * |Λ_d|^{1/2}
X_ase = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

% Procrustes alignment
[~, X_ase_aligned] = procrustes(X_true, X_ase);
sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);

fprintf('  ASE SSE: %.4f\n\n', sse_ase);

%% Surrogate Likelihood Optimization
fprintf('========================================\n');
fprintf('Surrogate Likelihood Optimization\n');
fprintf('========================================\n\n');

% Optimization options (matching max_iter = 30)
% TESTING: Using numerical gradients
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'iter', ...
    'MaxIterations', max_iter, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-8, ...
    'FunctionTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', false);  % Let MATLAB compute gradients

fprintf('Starting optimization...\n');
fprintf('  Initialization: ASE (Adjacency Spectral Embedding)\n');
fprintf('  X_init: %d x %d matrix from top %d eigenvectors\n', n, d, d);
fprintf('  Z_init: %d x %d matrix from least squares B ~ Z*X^T\n\n', p_cov, d);

% Add core folder to path
addpath('core');

tic;
[X_opt, Z_opt, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);
t_elapsed = toc;

%% Results
fprintf('\n========================================\n');
fprintf('Results\n');
fprintf('========================================\n\n');

% Procrustes alignment
[~, X_opt_aligned] = procrustes(X_true, X_opt);
sse_opt = sum((X_opt_aligned(:) - X_true(:)).^2);

fprintf('Convergence:\n');
fprintf('  Exit flag: %d\n', exitflag);
fprintf('  Iterations: %d / %d\n', output.iterations, max_iter);
fprintf('  Function evaluations: %d\n', output.funcCount);
fprintf('  First-order optimality: %.6e\n', output.firstorderopt);
fprintf('  Time elapsed: %.2f seconds\n\n', t_elapsed);

fprintf('Sum of Squared Errors (SSE):\n');
fprintf('  ASE:       %.4f\n', sse_ase);
fprintf('  Surrogate: %.4f\n', sse_opt);
fprintf('  Improvement: %.1f%%\n\n', 100 * (sse_ase - sse_opt) / sse_ase);

% Verify objective function
fprintf('Objective Function Components:\n');
fprintf('  Checking dimensions: B is %dx%d, X_opt is %dx%d, Z_opt is %dx%d\n', ...
    size(B,1), size(B,2), size(X_opt,1), size(X_opt,2), size(Z_opt,1), size(Z_opt,2));
Y_opt = X_opt * S;
S_mat = X_opt * Y_opt';
S_mat(1:n+1:end) = 0;
[psi_val, Psi_val, ~] = psi_functions(S_mat, tau);
net_obj = sum((A(:) - S_mat(:)) .* psi_val(:) + Psi_val(:));
B_pred = Z_opt * X_opt';
cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);
total_obj = net_obj + cov_obj;
fprintf('  Network term:   %+.6e\n', net_obj);
fprintf('  Covariate term: %+.6e\n', cov_obj);
fprintf('  Total (to max): %+.6e\n', total_obj);
fprintf('  fval (to min):  %+.6e\n', fval);
fprintf('  (Should match: fval = -total)\n\n');

%% Success criteria
fprintf('========================================\n');
if exitflag > 0
    fprintf('✓ Optimization converged\n');
else
    fprintf('⚠ Optimization did not converge (exitflag = %d)\n', exitflag);
end

if sse_opt < 0.5 * sse_ase
    fprintf('✓ Surrogate improved over ASE by >50%%\n');
else
    fprintf('⚠ Limited improvement (expected >50%%)\n');
end

if output.firstorderopt < 1e-3
    fprintf('✓ First-order optimality satisfied\n');
else
    fprintf('⚠ First-order optimality = %.6e (may need more iterations)\n', ...
        output.firstorderopt);
end

fprintf('========================================\n');
