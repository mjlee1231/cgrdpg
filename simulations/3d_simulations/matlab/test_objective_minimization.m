% Test objective function minimization
% Verify that optimization is actually minimizing (or equivalently maximizing the true objective)
clear; clc;

% Add core folder to path
addpath('core');

fprintf('Testing Objective Function Minimization\n');
fprintf('========================================\n\n');

%% Parameters
n = 100;
p_cov = 50;
d = 2;  % Changed to 2D
p = 2;  % All positive (no negative eigenvalues)
tau = 0.001;
rng(42);

fprintf('Parameters: n=%d, p_cov=%d, d=%d\n\n', n, p_cov, d);

%% Generate data (2D latent positions)
% theta <- pi * (1:n) / (n-1)
% X0[, 1] <- 0.28 * sin(theta) + 0.42
% X0[, 2] <- 0.28 * cos(theta) + 0.42
theta = pi * (1:n)' / (n - 1);
X_true = [0.28 * sin(theta) + 0.42, ...
          0.28 * cos(theta) + 0.42];

% Generate Probability Matrix P (True model)
P_net = X_true * X_true';

fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P_net(:)), max(P_net(:)));

% Generate Observed Adjacency Matrix A (Bernoulli)
A = double(rand(n) < P_net);
A = triu(A, 1);  % Upper triangle only
A = A + A';      % Make symmetric
A(1:n+1:end) = 0;  % Ensure diagonal is 0

% Generate covariates
Z_true = randn(p_cov, d) * 0.1;
B = Z_true * X_true' + randn(p_cov, n) * 0.1;

% For 2D case with all positive eigenvalues
S = eye(d);  % Identity signature (all positive)

%% Run optimization with iteration display
% Use 'iter' display to see objective values at each iteration
options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'iter', ...  % Show iteration details
    'MaxIterations', 50, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-8, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'lbfgs');

fprintf('Running fit_grdpg_fminunc_surrogate with Display=iter\n');
fprintf('Note: fval shown is NEGATIVE objective (for minimization)\n');
fprintf('The TRUE objective (to maximize) = -fval\n\n');

[X_opt, Z_opt, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);

%% Summary
fprintf('\n========================================\n');
fprintf('Optimization Summary\n');
fprintf('========================================\n\n');

fprintf('Exit flag: %d\n', exitflag);
fprintf('Total outer iterations: %d\n', output.iterations);
fprintf('Final fval (minimized): %+.6e\n', fval);
fprintf('Final objective (maximized): %+.6e\n', -fval);

fprintf('\nCheck the iteration output above:\n');
fprintf('  - Each outer iteration should show decreasing fval (for minimization)\n');
fprintf('  - Equivalently, -fval should be INCREASING (for maximization)\n');
fprintf('  - Inner iterations should show fval decreasing or staying similar\n');

%% Compute final objective components
Y_opt = X_opt * S;  % For 2D with identity S, Y_opt = X_opt
S_mat = X_opt * Y_opt';
is_off_diag = ~eye(n);

[psi_val, Psi_val, ~] = psi_functions(S_mat, tau);
net_obj = sum((A(is_off_diag) - S_mat(is_off_diag)) .* psi_val(is_off_diag) + Psi_val(is_off_diag));

B_pred = Z_opt * X_opt';
cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);

total_obj_max = net_obj + cov_obj;

fprintf('\nFinal objective breakdown:\n');
fprintf('  Network component:   %+.6e\n', net_obj);
fprintf('  Covariate component: %+.6e\n', cov_obj);
fprintf('  Total (to maximize): %+.6e\n', total_obj_max);
fprintf('  fval (to minimize):  %+.6e\n', -total_obj_max);
fprintf('  Reported fval:       %+.6e\n', fval);

if abs(fval - (-total_obj_max)) < 1e-6
    fprintf('  ✓ fval matches computed objective\n');
else
    fprintf('  ✗ WARNING: fval does not match!\n');
end
